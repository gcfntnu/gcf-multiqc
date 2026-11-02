"""Retain every index contribution in upstream BCL Convert sample/lane chunks."""

from collections import Counter, defaultdict
import csv

from multiqc.modules.bclconvert.bclconvert import ChunkMetrics, MultiqcModule as UpstreamBclConvert


def _join_indices(*values):
    return ", ".join(dict.fromkeys(index for value in values if value for index in value.split(", "))) or None


class MultiqcModule(UpstreamBclConvert):
    def parse_demux_data(self, demux_file, data_by_run, data_by_sample, num_demux_files):
        # Upstream totals are additive, but its chunk assignments overwrite earlier
        # index rows. Save prior chunks before it replaces them (also across files).
        previous = {}
        if demux_file.run_id in data_by_run:
            previous = {
                (lane_id, sample): chunk
                for lane_id, lane in data_by_run[demux_file.run_id].lanes.items()
                for sample, chunk in lane.samples.items()
            }
        super().parse_demux_data(demux_file, data_by_run, data_by_sample, num_demux_files)

        summed = {}
        with demux_file.path.open() as handle:
            for row in csv.DictReader(handle):
                sample = row["SampleID"]
                if sample == "Undetermined" or self.is_ignore_sample(sample):
                    continue
                key = (f"L{row['Lane']}", sample)
                if key not in summed:
                    summed[key] = (
                        previous[key].model_copy(deep=True)
                        if key in previous
                        else ChunkMetrics(run_id=demux_file.run_id, cluster_length=demux_file.cluster_length)
                    )
                chunk = summed[key]
                count = int(row["# Reads"])
                bases = count * demux_file.cluster_length
                chunk.clusters += count
                chunk.calculated_yield += bases
                chunk.perfect_index_reads += int(row["# Perfect Index Reads"])
                chunk.one_mismatch_index_reads += int(row["# One Mismatch Index Reads"])
                chunk.index = _join_indices(chunk.index, row["Index"])
                # Before 3.9.3 these fields lived in the demultiplexing CSV.
                if (q30 := row.get("# of >= Q30 Bases (PF)")) is not None:
                    chunk.yield_q30 = (chunk.yield_q30 or 0) + int(q30)
                if (quality := row.get("Mean Quality Score (PF)")) is not None:
                    chunk.calculated_qscore_sum = (chunk.calculated_qscore_sum or 0) + float(quality) * bases

        for (lane_id, sample), summed_chunk in summed.items():
            chunk = data_by_run[demux_file.run_id].lanes[lane_id].samples[sample]
            for field in (
                "clusters",
                "calculated_yield",
                "perfect_index_reads",
                "one_mismatch_index_reads",
                "yield_q30",
                "calculated_qscore_sum",
                "index",
            ):
                setattr(chunk, field, getattr(summed_chunk, field))

        # Share canonical chunks so later Quality_Metrics parsing updates both
        # views. Only colliding lane labels across runs need a run-name prefix.
        chunks_by_sample = defaultdict(list)
        for run_id, run in data_by_run.items():
            for lane_id, lane in run.lanes.items():
                for sample, chunk in lane.samples.items():
                    chunks_by_sample[sample].append((run_id, lane_id, chunk))
        for sample, chunks in chunks_by_sample.items():
            counts = Counter(lane_id for _, lane_id, _ in chunks)
            data_by_sample[sample].lanes = {
                f"{run_id} - {lane_id}" if counts[lane_id] > 1 else lane_id: chunk for run_id, lane_id, chunk in chunks
            }
            data_by_sample[sample].index = _join_indices(*(chunk.index for _, _, chunk in chunks))

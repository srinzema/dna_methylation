from pathlib import Path
import pandas as pd


def load_samples(sample_file: str, fastq_dir: str) -> pd.DataFrame:
    _df = pd.read_table(sample_file).set_index("samples", drop=True)
    if "alias" not in _df.columns:
        _df["alias"] = _df.index

    read1 = []
    read2 = []
    fastq_dir = Path(fastq_dir)
    for sample in _df.index:
        fastqs = [str(x) for x in fastq_dir.glob(f"{sample}*")]

        if len(fastqs) == 1:
            # one file found
            read1.append(fastqs[0])
            read2.append(pd.NA)
        elif len(fastqs) == 2:
            # two files found
            r1 = next(filter(lambda x: "R1" in x, fastqs))
            read1.append(r1)

            r2 = next(filter(lambda x: "R2" in x, fastqs))
            read2.append(r2)

    _df.insert(0, "read1", read1)
    _df.insert(1, "read2", read2)
    return _df

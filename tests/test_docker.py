from outward_assembly.io_helpers import s3_files_with_prefix
from outward_assembly.pipeline import outward_assembly
import os

paths = list(s3_files_with_prefix("lee-oa-test-data", "outward-assembly-test/siz/test_reads"))
print(f"Found {len(paths)} SIZ files")
for p in paths:
    print(f"  {p}")

work_dir = "./oa_workdir_docker"
os.makedirs(work_dir, exist_ok=True)

outward_assembly(
    seed_path="/workspace/tests/data/pipeline-end-to-end/seed.fasta",
    s3_paths=paths,
    output_path="./output_contigs.fasta",
    use_batch=False,
    work_dir_parent=work_dir,
    cleanup=False,
)

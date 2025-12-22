#!/usr/bin/env python3
"""
Generate SIZ-formatted test data for outward-assembly pipeline testing.

SIZ format:
- Split: chunked into files (default 1M read pairs, but we use smaller for testing)
- Interleaved: paired reads alternating <fwd1><rev1><fwd2><rev2>...
- Zstd-compressed: .fastq.zst extension

Naming convention: <prefix>_chunk000001.fastq.zst
"""

import argparse
import random
import subprocess
import tempfile
import os
from pathlib import Path


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))


def random_dna(length: int) -> str:
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choices('ATCG', k=length))


def generate_quality_string(length: int, min_qual: int = 30, max_qual: int = 40) -> str:
    """Generate a random quality string (Phred+33 encoded)."""
    return ''.join(chr(random.randint(min_qual, max_qual) + 33) for _ in range(length))


def generate_read_pair_with_seed(
    seed: str,
    read_length: int = 150,
    insert_size: int = 300,
    seed_position: str = "random"
) -> tuple[str, str, str, str]:
    """
    Generate a read pair where the forward read contains the seed.
    
    Returns: (fwd_seq, fwd_qual, rev_seq, rev_qual)
    """
    # Create a "fragment" that contains the seed
    # The fragment is insert_size long, with seed embedded somewhere
    
    if seed_position == "random":
        # Place seed at random position within the fragment
        max_pos = insert_size - len(seed)
        if max_pos <= 0:
            seed_start = 0
        else:
            seed_start = random.randint(0, max_pos)
    elif seed_position == "start":
        seed_start = 0
    elif seed_position == "middle":
        seed_start = (insert_size - len(seed)) // 2
    else:
        seed_start = 0
    
    # Build the fragment: random DNA + seed + random DNA
    prefix_len = seed_start
    suffix_len = insert_size - seed_start - len(seed)
    
    fragment = random_dna(prefix_len) + seed + random_dna(max(0, suffix_len))
    
    # Trim or pad to exact insert_size
    if len(fragment) < insert_size:
        fragment += random_dna(insert_size - len(fragment))
    fragment = fragment[:insert_size]
    
    # Forward read is first read_length bases
    fwd_seq = fragment[:read_length]
    fwd_qual = generate_quality_string(read_length)
    
    # Reverse read is reverse complement of last read_length bases
    rev_seq = reverse_complement(fragment[-read_length:])
    rev_qual = generate_quality_string(read_length)
    
    return fwd_seq, fwd_qual, rev_seq, rev_qual


def generate_random_read_pair(read_length: int = 150) -> tuple[str, str, str, str]:
    """Generate a completely random read pair (no seed)."""
    fwd_seq = random_dna(read_length)
    fwd_qual = generate_quality_string(read_length)
    rev_seq = random_dna(read_length)
    rev_qual = generate_quality_string(read_length)
    return fwd_seq, fwd_qual, rev_seq, rev_qual


def write_interleaved_fastq(
    filepath: Path,
    seed: str,
    num_pairs: int,
    seed_fraction: float = 0.1,
    read_length: int = 150,
    insert_size: int = 300
) -> int:
    """
    Write an interleaved FASTQ file with a mix of seed-containing and random reads.
    
    Returns: number of seed-containing read pairs written
    """
    seed_count = 0
    
    with open(filepath, 'w') as f:
        for i in range(num_pairs):
            read_id = f"read_{i+1:08d}"
            
            # Decide if this read pair should contain the seed
            if random.random() < seed_fraction:
                fwd_seq, fwd_qual, rev_seq, rev_qual = generate_read_pair_with_seed(
                    seed, read_length, insert_size
                )
                seed_count += 1
            else:
                fwd_seq, fwd_qual, rev_seq, rev_qual = generate_random_read_pair(read_length)
            
            # Write forward read
            f.write(f"@{read_id}/1\n")
            f.write(f"{fwd_seq}\n")
            f.write("+\n")
            f.write(f"{fwd_qual}\n")
            
            # Write reverse read (interleaved)
            f.write(f"@{read_id}/2\n")
            f.write(f"{rev_seq}\n")
            f.write("+\n")
            f.write(f"{rev_qual}\n")
    
    return seed_count


def compress_with_zstd(input_path: Path, output_path: Path, compression_level: int = 15) -> None:
    """Compress a file using zstd."""
    cmd = ["zstd", f"-{compression_level}", "-f", str(input_path), "-o", str(output_path)]
    subprocess.run(cmd, check=True)


def upload_to_s3(local_path: Path, s3_uri: str) -> None:
    """Upload a file to S3."""
    cmd = ["aws", "s3", "cp", str(local_path), s3_uri]
    subprocess.run(cmd, check=True)


def write_seed_fasta(filepath: Path, seed: str, seed_name: str = "seed") -> None:
    """Write seed sequence to a FASTA file."""
    with open(filepath, 'w') as f:
        f.write(f">{seed_name}\n")
        f.write(f"{seed}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate SIZ-formatted test data for outward-assembly"
    )
    parser.add_argument(
        "--seed",
        type=str,
        default="ATCGATCGATCGATCGATCGATCG",
        help="Seed sequence to embed in reads (default: ATCGATCGATCGATCGATCGATCG)"
    )
    parser.add_argument(
        "--num-pairs",
        type=int,
        default=10000,
        help="Number of read pairs per chunk (default: 10000)"
    )
    parser.add_argument(
        "--num-chunks",
        type=int,
        default=3,
        help="Number of SIZ chunks to generate (default: 3)"
    )
    parser.add_argument(
        "--seed-fraction",
        type=float,
        default=0.05,
        help="Fraction of reads that should contain the seed (default: 0.05)"
    )
    parser.add_argument(
        "--read-length",
        type=int,
        default=150,
        help="Read length in bp (default: 150)"
    )
    parser.add_argument(
        "--insert-size",
        type=int,
        default=300,
        help="Insert size in bp (default: 300)"
    )
    parser.add_argument(
        "--s3-bucket",
        type=str,
        required=True,
        help="S3 bucket name (without s3:// prefix)"
    )
    parser.add_argument(
        "--s3-prefix",
        type=str,
        default="outward-assembly-test/siz",
        help="S3 key prefix (default: outward-assembly-test/siz)"
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="test_reads",
        help="Output file prefix (default: test_reads)"
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=15,
        help="Zstd compression level (default: 15, matching NAO standards)"
    )
    parser.add_argument(
        "--local-output-dir",
        type=str,
        default=None,
        help="Local directory to keep output files (default: use temp dir and cleanup)"
    )
    parser.add_argument(
        "--skip-upload",
        action="store_true",
        help="Skip S3 upload (just generate local files)"
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="Random seed for reproducibility"
    )
    
    args = parser.parse_args()
    
    # Set random seed if provided
    if args.random_seed is not None:
        random.seed(args.random_seed)
    
    # Validate seed
    valid_bases = set('ATCGN')
    if not all(base.upper() in valid_bases for base in args.seed):
        raise ValueError(f"Seed contains invalid characters. Only ATCGN allowed.")
    
    seed = args.seed.upper()
    print(f"Seed sequence ({len(seed)} bp): {seed}")
    print(f"Generating {args.num_chunks} chunks with {args.num_pairs} read pairs each")
    print(f"Target seed fraction: {args.seed_fraction:.1%}")
    
    # Create output directory
    if args.local_output_dir:
        output_dir = Path(args.local_output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        output_dir = Path(tempfile.mkdtemp(prefix="siz_test_"))
        cleanup = True
    
    print(f"Working directory: {output_dir}")
    
    # Write seed FASTA file
    seed_fasta_path = output_dir / "seed.fasta"
    write_seed_fasta(seed_fasta_path, seed)
    print(f"Wrote seed FASTA: {seed_fasta_path}")
    
    # Generate chunks
    total_seed_reads = 0
    generated_files = []
    
    for chunk_idx in range(1, args.num_chunks + 1):
        chunk_name = f"{args.output_prefix}_chunk{chunk_idx:06d}"
        fastq_path = output_dir / f"{chunk_name}.fastq"
        zst_path = output_dir / f"{chunk_name}.fastq.zst"
        
        print(f"\nGenerating chunk {chunk_idx}/{args.num_chunks}: {chunk_name}")
        
        # Generate FASTQ
        seed_count = write_interleaved_fastq(
            fastq_path,
            seed=seed,
            num_pairs=args.num_pairs,
            seed_fraction=args.seed_fraction,
            read_length=args.read_length,
            insert_size=args.insert_size
        )
        total_seed_reads += seed_count
        print(f"  - Wrote {args.num_pairs} read pairs ({seed_count} contain seed)")
        
        # Compress with zstd
        compress_with_zstd(fastq_path, zst_path, args.compression_level)
        fastq_size = fastq_path.stat().st_size
        zst_size = zst_path.stat().st_size
        print(f"  - Compressed: {fastq_size:,} -> {zst_size:,} bytes ({zst_size/fastq_size:.1%})")
        
        # Remove uncompressed file to save space
        fastq_path.unlink()
        
        generated_files.append(zst_path)
    
    print(f"\n{'='*60}")
    print(f"Generated {args.num_chunks} SIZ files")
    print(f"Total read pairs: {args.num_pairs * args.num_chunks:,}")
    print(f"Total seed-containing pairs: {total_seed_reads:,} ({total_seed_reads/(args.num_pairs * args.num_chunks):.2%})")
    
    # Upload to S3
    if not args.skip_upload:
        print(f"\nUploading to s3://{args.s3_bucket}/{args.s3_prefix}/")
        
        # Upload seed FASTA
        seed_s3_uri = f"s3://{args.s3_bucket}/{args.s3_prefix}/seed.fasta"
        upload_to_s3(seed_fasta_path, seed_s3_uri)
        print(f"  - Uploaded seed: {seed_s3_uri}")
        
        # Upload SIZ files
        s3_paths = []
        for zst_path in generated_files:
            s3_uri = f"s3://{args.s3_bucket}/{args.s3_prefix}/{zst_path.name}"
            upload_to_s3(zst_path, s3_uri)
            s3_paths.append(s3_uri)
            print(f"  - Uploaded: {s3_uri}")
        
        print(f"\n{'='*60}")
        print("S3 PATHS FOR OUTWARD ASSEMBLY:")
        print(f"{'='*60}")
        print(f"Seed FASTA: {seed_s3_uri}")
        print(f"\nSIZ files prefix: s3://{args.s3_bucket}/{args.s3_prefix}/{args.output_prefix}")
        print("\nTo use with outward_assembly:")
        print(f'''
from outward_assembly.io_helpers import s3_files_with_prefix
from outward_assembly.pipeline import outward_assembly

# Get S3 paths
paths = s3_files_with_prefix("{args.s3_bucket}", "{args.s3_prefix}/{args.output_prefix}")

# Run assembly (local profile)
outward_assembly(
    seed_path="seed.fasta",  # Download from S3 or use local copy
    s3_paths=paths,
    output_path="output_contigs.fasta",
    use_batch=False,
)
''')
    else:
        print("\nSkipped S3 upload (--skip-upload flag set)")
        print(f"Local files in: {output_dir}")
    
    # Cleanup temp directory if needed
    if cleanup and not args.skip_upload:
        import shutil
        shutil.rmtree(output_dir)
        print(f"\nCleaned up temporary directory")
    elif not cleanup:
        print(f"\nFiles retained in: {output_dir}")


if __name__ == "__main__":
    main()
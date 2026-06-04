#!/usr/bin/env python3
"""Find pneumococcal cps genes in assemblies, extract hit sequences, and align by gene.

The workflow is intentionally lightweight: it uses BLAST+ and MAFFT externally and
keeps parsing/extraction in the Python standard library so it can run on clusters
without extra Python packages.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

FASTA_EXTS = (".fa", ".fna", ".fasta", ".fas", ".fa.gz", ".fna.gz", ".fasta.gz")
BLAST_OUTFMT = "6 qseqid sseqid qstart qend sstart send pident length evalue bitscore qlen slen"
DEFAULT_GAP_N = 20


@dataclass(frozen=True)
class FastaRecord:
    name: str
    description: str
    sequence: str


@dataclass
class Hit:
    isolate: str
    gene: str
    source: str
    contig: str
    qstart: int
    qend: int
    sstart: int
    send: int
    pident: float
    length: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int

    @property
    def qlo(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def qhi(self) -> int:
        return max(self.qstart, self.qend)

    @property
    def slo(self) -> int:
        return min(self.sstart, self.send)

    @property
    def shi(self) -> int:
        return max(self.sstart, self.send)

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"


@dataclass
class GeneCall:
    isolate: str
    gene: str
    source: str
    status: str
    qlen: int
    qcov: float
    mean_identity: float
    total_bitscore: float
    n_hsps: int
    n_contigs: int
    max_query_gap: int
    notes: str
    hits: list[Hit]


@dataclass(frozen=True)
class BreakEvent:
    isolate: str
    gene: str
    status: str
    source: str
    break_type: str
    query_left_end: int
    query_right_start: int
    missing_query_length: int
    left_hsp: str
    right_hsp: str
    left_contig: str
    left_subject_start: str
    left_subject_end: str
    right_contig: str
    right_subject_start: str
    right_subject_end: str
    strands: str
    assembly_gap_bases: str
    assembly_gap_delta_from_query: str
    assembly_gap_mod3: str
    frameshift_suspect: str
    notes: str


def log(message: str) -> None:
    print(f"[cps-pipeline] {message}", file=sys.stderr, flush=True)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def read_fasta(path: Path) -> Iterator[FastaRecord]:
    name = None
    desc = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield FastaRecord(name, desc, "".join(chunks).upper())
                desc = line[1:]
                name = desc.split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield FastaRecord(name, desc, "".join(chunks).upper())


def write_fasta(records: Iterable[tuple[str, str]], path: Path, width: int = 80) -> None:
    with path.open("w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
    return seq.translate(table)[::-1].upper()


def discover_fastas(directory: Path) -> list[Path]:
    files = [p for p in directory.iterdir() if p.is_file() and str(p).lower().endswith(FASTA_EXTS)]
    return sorted(files)


def require_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"Required executable not found on PATH: {name}")


def run_command(cmd: list[str], dry_run: bool = False) -> None:
    log(" ".join(cmd))
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def make_blast_db(assembly: Path, db_prefix: Path, dry_run: bool = False) -> None:
    nin = db_prefix.with_suffix(".nin")
    nsq = db_prefix.with_suffix(".nsq")
    if nin.exists() or nsq.exists():
        return
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    run_command(["makeblastdb", "-in", str(assembly), "-dbtype", "nucl", "-out", str(db_prefix)], dry_run)


def blast_query(tool: str, query: Path, db_prefix: Path, out_path: Path, threads: int, evalue: str, dry_run: bool = False) -> None:
    if out_path.exists() and out_path.stat().st_size > 0:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        tool,
        "-query",
        str(query),
        "-db",
        str(db_prefix),
        "-outfmt",
        BLAST_OUTFMT,
        "-evalue",
        evalue,
        "-num_threads",
        str(threads),
        "-out",
        str(out_path),
    ]
    if tool == "blastn":
        cmd.extend(["-task", "blastn"])
    run_command(cmd, dry_run)


def parse_blast(path: Path, isolate: str, source: str, min_identity: float, min_hsp_len: int) -> Iterator[Hit]:
    if not path.exists():
        return
    with path.open() as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 12:
                raise ValueError(f"Expected 12 BLAST columns in {path}, saw {len(fields)}: {raw[:120]}")
            qseqid, sseqid, qstart, qend, sstart, send, pident, length, evalue, bitscore, qlen, slen = fields
            hit = Hit(
                isolate=isolate,
                gene=qseqid,
                source=source,
                contig=sseqid,
                qstart=int(qstart),
                qend=int(qend),
                sstart=int(sstart),
                send=int(send),
                pident=float(pident),
                length=int(length),
                evalue=float(evalue),
                bitscore=float(bitscore),
                qlen=int(qlen),
                slen=int(slen),
            )
            if hit.pident >= min_identity and hit.length >= min_hsp_len:
                yield hit


def greedy_non_overlapping_hits(hits: list[Hit], max_query_overlap: int) -> list[Hit]:
    """Keep a best-scoring set of HSPs that together cover the query."""
    chosen: list[Hit] = []
    for hit in sorted(hits, key=lambda h: (-h.bitscore, -(h.qhi - h.qlo + 1), h.qlo)):
        overlaps_too_much = False
        for existing in chosen:
            overlap = max(0, min(hit.qhi, existing.qhi) - max(hit.qlo, existing.qlo) + 1)
            if overlap > max_query_overlap:
                overlaps_too_much = True
                break
        if not overlaps_too_much:
            chosen.append(hit)
    return sorted(chosen, key=lambda h: (h.qlo, h.qhi, h.contig, h.slo))


def covered_query_bases(hits: list[Hit]) -> int:
    intervals = sorted((h.qlo, h.qhi) for h in hits)
    if not intervals:
        return 0
    merged: list[list[int]] = []
    for lo, hi in intervals:
        if not merged or lo > merged[-1][1] + 1:
            merged.append([lo, hi])
        else:
            merged[-1][1] = max(merged[-1][1], hi)
    return sum(hi - lo + 1 for lo, hi in merged)


def max_query_gap(hits: list[Hit]) -> int:
    ordered = sorted(hits, key=lambda h: h.qlo)
    if len(ordered) < 2:
        return 0
    return max(max(0, nxt.qlo - prev.qhi - 1) for prev, nxt in zip(ordered, ordered[1:]))


def classify_hits(
    hits: list[Hit],
    complete_qcov: float,
    partial_qcov: float,
    fragmented_gap: int,
    max_query_overlap: int,
) -> GeneCall | None:
    if not hits:
        return None

    selected = greedy_non_overlapping_hits(hits, max_query_overlap=max_query_overlap)
    qlen = max(h.qlen for h in selected)
    qcov = covered_query_bases(selected) / qlen if qlen else 0.0
    weighted_identity = sum(h.pident * h.length for h in selected) / max(1, sum(h.length for h in selected))
    n_contigs = len({h.contig for h in selected})
    qgap = max_query_gap(selected)
    fragmented = len(selected) > 1 and (n_contigs > 1 or qgap >= fragmented_gap)

    if qcov >= complete_qcov and not fragmented:
        status = "complete"
        notes = "single high-coverage hit"
    elif qcov >= complete_qcov and fragmented:
        status = "fragmented"
        notes = "high coverage split across HSPs/contigs"
    elif qcov >= partial_qcov:
        status = "partial"
        notes = "sub-threshold query coverage"
    else:
        status = "absent"
        notes = "only low-coverage hits found"

    return GeneCall(
        isolate=selected[0].isolate,
        gene=selected[0].gene,
        source=selected[0].source,
        status=status,
        qlen=qlen,
        qcov=qcov,
        mean_identity=weighted_identity,
        total_bitscore=sum(h.bitscore for h in selected),
        n_hsps=len(selected),
        n_contigs=n_contigs,
        max_query_gap=qgap,
        notes=notes,
        hits=selected,
    )


def best_call(calls: list[GeneCall]) -> GeneCall | None:
    if not calls:
        return None
    status_rank = {"complete": 4, "fragmented": 3, "partial": 2, "absent": 1}
    return max(calls, key=lambda c: (status_rank[c.status], c.qcov, c.total_bitscore, c.mean_identity))


def load_assembly(path: Path) -> dict[str, str]:
    return {record.name: record.sequence for record in read_fasta(path)}


def extract_call_sequence(call: GeneCall, contigs: dict[str, str], gap_n: int) -> str:
    pieces: list[str] = []
    for hit in sorted(call.hits, key=lambda h: (h.qlo, h.qhi)):
        try:
            seq = contigs[hit.contig][hit.slo - 1 : hit.shi]
        except KeyError as exc:
            raise KeyError(f"Contig {hit.contig!r} from BLAST output not found in assembly for {call.isolate}") from exc
        if hit.strand == "-":
            seq = revcomp(seq)
        pieces.append(seq.upper())
    return ("N" * gap_n).join(pieces)


def expected_subject_gap(source: str, query_gap_length: int) -> int:
    if source == "tblastn":
        return query_gap_length * 3
    return query_gap_length


def subject_gap_between(left: Hit, right: Hit) -> int | None:
    if left.contig != right.contig or left.strand != right.strand:
        return None
    if left.strand == "+":
        return right.slo - left.shi - 1
    return left.slo - right.shi - 1


def break_type_for_hsp_pair(left: Hit, right: Hit, query_gap_length: int, fragmented_gap: int, frameshift_suspect: bool) -> str:
    if left.contig != right.contig:
        return "contig_break"
    if left.strand != right.strand:
        return "strand_switch"
    if frameshift_suspect:
        return "possible_frameshift"
    if query_gap_length >= fragmented_gap:
        return "internal_query_gap"
    return "hsp_split"


def break_events_for_call(call: GeneCall, fragmented_gap: int) -> list[BreakEvent]:
    ordered = sorted(call.hits, key=lambda h: (h.qlo, h.qhi))
    if not ordered:
        return []

    events: list[BreakEvent] = []
    first = ordered[0]
    last = ordered[-1]

    if first.qlo > 1:
        events.append(
            BreakEvent(
                isolate=call.isolate,
                gene=call.gene,
                status=call.status,
                source=call.source,
                break_type="missing_5_prime",
                query_left_end=0,
                query_right_start=first.qlo,
                missing_query_length=first.qlo - 1,
                left_hsp="",
                right_hsp="1",
                left_contig="",
                left_subject_start="",
                left_subject_end="",
                right_contig=first.contig,
                right_subject_start=str(first.sstart),
                right_subject_end=str(first.send),
                strands=first.strand,
                assembly_gap_bases="",
                assembly_gap_delta_from_query="",
                assembly_gap_mod3="",
                frameshift_suspect="no",
                notes="reference/query start is not covered by selected HSPs",
            )
        )

    for idx, (left, right) in enumerate(zip(ordered, ordered[1:]), start=1):
        query_gap_length = max(0, right.qlo - left.qhi - 1)
        assembly_gap = subject_gap_between(left, right)
        expected_gap = expected_subject_gap(call.source, query_gap_length)
        delta = assembly_gap - expected_gap if assembly_gap is not None else None
        frameshift_suspect = call.source == "tblastn" and delta is not None and delta % 3 != 0
        break_type = break_type_for_hsp_pair(left, right, query_gap_length, fragmented_gap, frameshift_suspect)
        if assembly_gap is None:
            gap_text = ""
            delta_text = ""
            mod_text = ""
        else:
            gap_text = str(assembly_gap)
            delta_text = str(delta)
            mod_text = str(abs(delta) % 3)

        notes = []
        if assembly_gap is not None and assembly_gap < 0:
            notes.append("selected HSPs overlap on the assembly")
        if assembly_gap is not None and assembly_gap > 0:
            notes.append("unmatched assembly sequence lies between selected HSPs")
        if query_gap_length > 0:
            notes.append("reference/query positions are missing between selected HSPs")
        if frameshift_suspect:
            notes.append("tBLASTn split has a non-triplet assembly/query gap delta")
        if not notes:
            notes.append("selected HSPs are adjacent or nearly adjacent")

        events.append(
            BreakEvent(
                isolate=call.isolate,
                gene=call.gene,
                status=call.status,
                source=call.source,
                break_type=break_type,
                query_left_end=left.qhi,
                query_right_start=right.qlo,
                missing_query_length=query_gap_length,
                left_hsp=str(idx),
                right_hsp=str(idx + 1),
                left_contig=left.contig,
                left_subject_start=str(left.sstart),
                left_subject_end=str(left.send),
                right_contig=right.contig,
                right_subject_start=str(right.sstart),
                right_subject_end=str(right.send),
                strands=f"{left.strand}/{right.strand}",
                assembly_gap_bases=gap_text,
                assembly_gap_delta_from_query=delta_text,
                assembly_gap_mod3=mod_text,
                frameshift_suspect="yes" if frameshift_suspect else "no",
                notes="; ".join(notes),
            )
        )

    if last.qhi < call.qlen:
        events.append(
            BreakEvent(
                isolate=call.isolate,
                gene=call.gene,
                status=call.status,
                source=call.source,
                break_type="missing_3_prime",
                query_left_end=last.qhi,
                query_right_start=call.qlen + 1,
                missing_query_length=call.qlen - last.qhi,
                left_hsp=str(len(ordered)),
                right_hsp="",
                left_contig=last.contig,
                left_subject_start=str(last.sstart),
                left_subject_end=str(last.send),
                right_contig="",
                right_subject_start="",
                right_subject_end="",
                strands=last.strand,
                assembly_gap_bases="",
                assembly_gap_delta_from_query="",
                assembly_gap_mod3="",
                frameshift_suspect="no",
                notes="reference/query end is not covered by selected HSPs",
            )
        )

    return events


def write_summary(calls: list[GeneCall], genes: list[str], isolates: list[str], path: Path) -> None:
    by_key = {(c.isolate, c.gene): c for c in calls}
    fields = [
        "isolate",
        "gene",
        "status",
        "source",
        "query_length",
        "query_coverage",
        "mean_identity",
        "total_bitscore",
        "hsp_count",
        "contig_count",
        "max_query_gap",
        "notes",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for isolate in isolates:
            for gene in genes:
                call = by_key.get((isolate, gene))
                if call is None:
                    row = {
                        "isolate": isolate,
                        "gene": gene,
                        "status": "absent",
                        "source": "none",
                        "query_length": "",
                        "query_coverage": "0.0000",
                        "mean_identity": "",
                        "total_bitscore": "0.0",
                        "hsp_count": "0",
                        "contig_count": "0",
                        "max_query_gap": "",
                        "notes": "no passing BLAST hits",
                    }
                else:
                    row = {
                        "isolate": call.isolate,
                        "gene": call.gene,
                        "status": call.status,
                        "source": call.source,
                        "query_length": str(call.qlen),
                        "query_coverage": f"{call.qcov:.4f}",
                        "mean_identity": f"{call.mean_identity:.2f}",
                        "total_bitscore": f"{call.total_bitscore:.1f}",
                        "hsp_count": str(call.n_hsps),
                        "contig_count": str(call.n_contigs),
                        "max_query_gap": str(call.max_query_gap),
                        "notes": call.notes,
                    }
                writer.writerow(row)


def write_bed(calls: list[GeneCall], path: Path) -> None:
    with path.open("w") as handle:
        for call in sorted(calls, key=lambda c: (c.isolate, c.gene)):
            for idx, hit in enumerate(call.hits, start=1):
                # BED is 0-based, half-open. BLAST coordinates are 1-based inclusive.
                start0 = hit.slo - 1
                end0 = hit.shi
                name = f"{call.isolate}|{call.gene}|{call.status}|{call.source}|hsp{idx}"
                handle.write(
                    "\t".join(
                        [
                            hit.contig,
                            str(start0),
                            str(end0),
                            name,
                            f"{hit.bitscore:.1f}",
                            hit.strand,
                        ]
                    )
                    + "\n"
                )


def write_breakpoints(calls: list[GeneCall], path: Path, fragmented_gap: int) -> None:
    fields = list(BreakEvent.__dataclass_fields__)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for call in sorted(calls, key=lambda c: (c.isolate, c.gene)):
            for event in break_events_for_call(call, fragmented_gap):
                writer.writerow({field: getattr(event, field) for field in fields})


def run_alignments(per_gene_dir: Path, aln_dir: Path, threads: int, dry_run: bool = False) -> None:
    require_executable("mafft")
    aln_dir.mkdir(parents=True, exist_ok=True)
    for fasta in sorted(per_gene_dir.glob("*.fna")):
        out_path = aln_dir / f"{fasta.stem}.aln.fna"
        cmd = ["mafft", "--thread", str(threads), "--auto", str(fasta)]
        log(" ".join(cmd) + f" > {out_path}")
        if dry_run:
            continue
        with out_path.open("w") as out_handle:
            subprocess.run(cmd, check=True, stdout=out_handle)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract and align cps genes from bacterial assemblies using BLAST+ and MAFFT.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--assemblies", type=Path, default=Path("inputs/assemblies"), help="Directory of isolate assemblies")
    parser.add_argument("--nucl-refs", type=Path, default=Path("inputs/refs/cps_genes.fna"), help="Nucleotide gene reference FASTA")
    parser.add_argument("--prot-refs", type=Path, default=Path("inputs/refs/cps_genes.faa"), help="Protein gene reference FASTA")
    parser.add_argument("--outdir", type=Path, default=Path("results/cps"), help="Output directory")
    parser.add_argument("--mode", choices=("blastn", "tblastn", "both"), default="both", help="Search strategy")
    parser.add_argument("--threads", type=int, default=1, help="Threads passed to BLAST and MAFFT")
    parser.add_argument("--evalue", default="1e-20", help="BLAST e-value threshold")
    parser.add_argument("--min-identity", type=float, default=70.0, help="Minimum HSP percent identity to keep")
    parser.add_argument("--min-hsp-len", type=int, default=30, help="Minimum HSP length to keep (nt for blastn, aa for tblastn)")
    parser.add_argument("--complete-qcov", type=float, default=0.90, help="Query coverage threshold for complete calls")
    parser.add_argument("--partial-qcov", type=float, default=0.20, help="Query coverage threshold for partial calls")
    parser.add_argument("--fragmented-gap", type=int, default=30, help="Query gap size that marks a multi-HSP hit as fragmented")
    parser.add_argument("--max-query-overlap", type=int, default=10, help="Allowed overlap when choosing non-overlapping HSPs")
    parser.add_argument("--gap-n", type=int, default=DEFAULT_GAP_N, help="Number of Ns inserted between fragmented HSPs")
    parser.add_argument("--skip-align", action="store_true", help="Do not run MAFFT")
    parser.add_argument("--dry-run", action="store_true", help="Print external commands without executing them")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    require_executable("makeblastdb")
    if args.mode in {"blastn", "both"}:
        require_executable("blastn")
    if args.mode in {"tblastn", "both"}:
        require_executable("tblastn")

    assemblies = discover_fastas(args.assemblies)
    if not assemblies:
        raise SystemExit(f"No assemblies found in {args.assemblies}")

    nucl_genes: list[str] = []
    prot_genes: list[str] = []
    if args.mode in {"blastn", "both"}:
        if not args.nucl_refs.exists():
            raise SystemExit(f"Nucleotide references not found: {args.nucl_refs}")
        nucl_genes = [r.name for r in read_fasta(args.nucl_refs)]
    if args.mode in {"tblastn", "both"}:
        if not args.prot_refs.exists():
            raise SystemExit(f"Protein references not found: {args.prot_refs}")
        prot_genes = [r.name for r in read_fasta(args.prot_refs)]
    genes = sorted(set(nucl_genes) | set(prot_genes))

    db_dir = args.outdir / "blastdb"
    blast_dir = args.outdir / "blast"
    seq_dir = args.outdir / "sequences"
    per_gene_dir = args.outdir / "per_gene"
    aln_dir = args.outdir / "alignments"
    for directory in (db_dir, blast_dir, seq_dir, per_gene_dir):
        directory.mkdir(parents=True, exist_ok=True)

    hits_by_key_source: dict[tuple[str, str, str], list[Hit]] = defaultdict(list)
    isolate_names: list[str] = []
    assembly_by_isolate: dict[str, Path] = {}

    for assembly in assemblies:
        isolate = assembly.name
        for ext in FASTA_EXTS:
            if isolate.lower().endswith(ext):
                isolate = isolate[: -len(ext)]
                break
        isolate_names.append(isolate)
        assembly_by_isolate[isolate] = assembly
        db_prefix = db_dir / isolate / "assembly"
        make_blast_db(assembly, db_prefix, args.dry_run)
        if args.mode in {"blastn", "both"}:
            out_path = blast_dir / isolate / "blastn.tsv"
            blast_query("blastn", args.nucl_refs, db_prefix, out_path, args.threads, args.evalue, args.dry_run)
            for hit in parse_blast(out_path, isolate, "blastn", args.min_identity, args.min_hsp_len):
                hits_by_key_source[(isolate, hit.gene, "blastn")].append(hit)
        if args.mode in {"tblastn", "both"}:
            out_path = blast_dir / isolate / "tblastn.tsv"
            blast_query("tblastn", args.prot_refs, db_prefix, out_path, args.threads, args.evalue, args.dry_run)
            for hit in parse_blast(out_path, isolate, "tblastn", args.min_identity, args.min_hsp_len):
                hits_by_key_source[(isolate, hit.gene, "tblastn")].append(hit)

    calls: list[GeneCall] = []
    for isolate in isolate_names:
        for gene in genes:
            candidate_calls = []
            for source in ("blastn", "tblastn"):
                call = classify_hits(
                    hits_by_key_source.get((isolate, gene, source), []),
                    complete_qcov=args.complete_qcov,
                    partial_qcov=args.partial_qcov,
                    fragmented_gap=args.fragmented_gap,
                    max_query_overlap=args.max_query_overlap,
                )
                if call is not None:
                    candidate_calls.append(call)
            call = best_call(candidate_calls)
            if call is not None and call.status != "absent":
                calls.append(call)

    extracted_records_by_gene: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for isolate in isolate_names:
        contigs = load_assembly(assembly_by_isolate[isolate])
        isolate_calls = [c for c in calls if c.isolate == isolate]
        isolate_records: list[tuple[str, str]] = []
        for call in sorted(isolate_calls, key=lambda c: c.gene):
            seq = extract_call_sequence(call, contigs, args.gap_n)
            header = f"{call.isolate}|{call.gene}|{call.status}|{call.source}|qcov={call.qcov:.3f}|pid={call.mean_identity:.1f}"
            isolate_records.append((header, seq))
            extracted_records_by_gene[call.gene].append((header, seq))
        write_fasta(isolate_records, seq_dir / f"{isolate}.cps_hits.fna")

    for gene, records in sorted(extracted_records_by_gene.items()):
        safe_gene = gene.replace(os.sep, "_")
        write_fasta(records, per_gene_dir / f"{safe_gene}.fna")

    write_summary(calls, genes, isolate_names, args.outdir / "summary.tsv")
    write_bed(calls, args.outdir / "cps_hits.bed")
    write_breakpoints(calls, args.outdir / "breakpoints.tsv", args.fragmented_gap)

    if not args.skip_align:
        run_alignments(per_gene_dir, aln_dir, args.threads, args.dry_run)

    log(f"Wrote outputs under {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

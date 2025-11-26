"""CLI for vcf2tsv - Convert VCF files to TSV format."""

import os
import re
import sys
from enum import Enum
from signal import SIG_DFL, SIGPIPE, signal
from subprocess import PIPE, Popen

import typer
from cyvcf2 import VCF
from typing_extensions import Annotated

# Handle broken pipe gracefully
try:
    signal(SIGPIPE, SIG_DFL)
except (AttributeError, ValueError):
    # SIGPIPE is not available on Windows
    pass

app = typer.Typer(
    name="vcf2tsv",
    help="Convert VCF files to TSV format.",
    add_completion=False,
)


class OutputFormat(str, Enum):
    wide = "wide"
    long = "long"


# Regex patterns for parsing VCF header
RE_INFO = re.compile(
    r"""^\#\#INFO=<
    ID=(?P<id>[^,]+),
    Number=(?P<number>-?\d+|\.|[AG]),
    Type=(?P<type>Integer|Float|Flag|Character|String),
    Description="(?P<desc>[^"]*)".*
    >""",
    re.VERBOSE | re.MULTILINE,
)

RE_FORMAT = re.compile(
    r"""^\#\#FORMAT=<
    ID=(?P<id>.+),
    Number=(?P<number>-?\d+|\.|[AGR]),
    Type=(?P<type>.+),
    Description="(?P<desc>.*)".*
    >""",
    re.VERBOSE | re.MULTILINE,
)

ANN_HEADER = [
    "allele",
    "effect",
    "impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "exon_intron_rank",
    "nt_change",
    "aa_change",
    "cDNA_position/cDNA_len",
    "protein_position",
    "distance_to_feature",
    "error",
]


def convert_vcf_to_tsv(
    vcf_path: str,
    output_format: OutputFormat,
    print_header: bool = False,
    expand_ann: bool = False,
) -> None:
    """Convert VCF to TSV format and print to stdout.

    Args:
        vcf_path: Path to the VCF file.
        output_format: Output format (wide or long).
        print_header: Whether to print the header line.
        expand_ann: Whether to expand snpEff ANN annotations.
    """
    if not os.path.isfile(vcf_path) and vcf_path != "-":
        typer.echo(f"Error: {vcf_path} does not exist", err=True)
        raise typer.Exit(1)

    vcf = VCF(vcf_path)
    raw_header = vcf.raw_header
    samples = vcf.samples
    filename = vcf_path

    # Parse INFO and FORMAT fields from header
    info_fields = [m.groupdict()["id"] for m in RE_INFO.finditer(raw_header)]
    format_fields = [m.groupdict()["id"] for m in RE_FORMAT.finditer(raw_header)]

    # Construct bcftools query string
    header_flag = "--print-header" if print_header else ""

    if output_format == OutputFormat.long:
        query_string = (
            "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t"
            + "\\t".join([f"%{x}" for x in info_fields])
            + "\\t[-->%SAMPLE\\t"
            + "\\t".join([f"%{x}" for x in format_fields])
            + "\\n]"
        )
    else:  # wide
        query_string = (
            "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t"
            + "\\t".join([f"%INFO/{x}" for x in info_fields])
            + "[\\t%SAMPLE\\t"
            + "\\t".join([f"%{x}" for x in format_fields])
            + "]\\n"
        )

    # Build bcftools command
    cmd = list(filter(len, ["bcftools", "query", header_flag, "-f", query_string, filename]))
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE)

    # Calculate ANN field position if needed
    ann_loc = None
    if expand_ann and "ANN" in info_fields:
        ann_loc = info_fields.index("ANN") + 7  # 7 = CHROM, POS, ID, REF, ALT, QUAL, FILTER

    fill_fields = []

    for n, line in enumerate(proc.stdout):
        line = line.decode("UTF-8")

        if n == 0 and print_header and output_format == OutputFormat.wide:
            # Wide format header
            if ann_loc is not None:
                parts = line.split("\t")
                parts = parts[: ann_loc - 1] + ANN_HEADER + parts[ann_loc + 1 :]
                line = "\t".join(parts)
            print(re.sub(r"\[[0-9]+\]", "", line).strip("#\n ").replace(":", "_"))

        elif n == 0 and print_header and output_format == OutputFormat.long:
            # Long format header
            line = re.sub(r"\[[0-9]+\]", "", line).strip("#\n ").split("\t")
            header = []
            for var in line:
                if ":" in var:
                    var = var.split(":")[1]
                header.append(var)
            if ann_loc is not None:
                header = header[: ann_loc - 1] + ANN_HEADER + header[ann_loc + 1 :]
            # Prefix FORMAT fields with F_
            header = header[: -len(format_fields)] + [f"F_{x}" for x in header[-len(format_fields) :]]
            print("\t".join(header).replace("-->", ""))

        elif n < len(samples) and output_format == OutputFormat.long:
            # Skip initial sample lines in long format
            pass

        elif n >= len(samples) and output_format == OutputFormat.long:
            # Long format data lines
            parts = line.strip("\n").split("\t")
            if parts[0].startswith("-->"):
                parts = fill_fields + parts[len(parts) - len(fill_fields) :]
            else:
                fill_fields = parts[0 : (7 + len(info_fields))]

            if ann_loc is not None:
                for var_effect in parts[ann_loc].split(","):
                    out_line = parts[: ann_loc - 1] + var_effect.split("|") + parts[ann_loc + 2 :]
                    print("\t".join(out_line).strip("\n").replace("-->", ""))
            else:
                print("\t".join(parts).replace("-->", ""))

        else:
            # Wide format data lines
            if ann_loc is not None:
                parts = line.split("\t")
                for var_effect in parts[ann_loc].split(","):
                    out_line = parts[: ann_loc - 1] + var_effect.split("|") + parts[ann_loc + 2 :]
                    print("\t".join(out_line).strip("\n"))
            else:
                print(line.strip("\n"))


@app.command()
def main(
    vcf: Annotated[str, typer.Argument(help="Path to VCF file")],
    format: Annotated[
        OutputFormat,
        typer.Option("--format", "-f", help="Output format"),
    ] = OutputFormat.wide,
    header: Annotated[
        bool,
        typer.Option("--header", "-h", help="Print header line"),
    ] = False,
    ann: Annotated[
        bool,
        typer.Option("--ann", "-a", help="Expand snpEff ANN annotations"),
    ] = False,
) -> None:
    """Convert VCF file to TSV format.

    Requires bcftools to be installed and available in PATH.
    """
    convert_vcf_to_tsv(vcf, format, header, ann)


if __name__ == "__main__":
    app()

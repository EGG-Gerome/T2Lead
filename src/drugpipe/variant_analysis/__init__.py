"""Somatic variant analysis: tumor/normal paired variant calling bridge.

Connects upstream genomic variant detection (BWA-MEM2 → Mutect2 → VEP)
to downstream structure-based drug discovery (ESMFold → Docking → MD).
"""

from drugpipe.variant_analysis.vcf_parser import VCFParser
from drugpipe.variant_analysis.mutant_sequence import MutantSequenceBuilder
from drugpipe.variant_analysis.structure_bridge import StructureBridge

__all__ = ["VCFParser", "MutantSequenceBuilder", "StructureBridge"]

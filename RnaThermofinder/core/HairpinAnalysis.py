from pathlib import Path

import RNA
from typing import List, Tuple, Callable, Optional, Dict, Any
import csv

# ✨ NEW: Import for composition and CSV building
import sys

from RnaThermofinder.utils.analysis_helpers import calculate_composition
from RnaThermofinder.utils.analysis_helpers import build_csv_row
from settings_manager import SettingsManager

sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def get_terminal_hairpin_with_tail(sequence, structure):
    """
    Extracts the rightmost (terminal) hairpin and all trailing unpaired dots.
    Example: ((((((((((...)))))))))))..........
    """
    # Find the rightmost ')'
    last_close = structure.rfind(')')
    if last_close == -1:
        return None  # no paired region found

    # Now find its matching '(' going backward
    depth = 0
    for i in range(last_close, -1, -1):
        if structure[i] == ')':
            depth += 1
        elif structure[i] == '(':
            depth -= 1
            if depth == 0:
                start = i
                break
    else:
        return None  # no matching '(' found

    # Include all trailing unpaired dots after the hairpin
    end = last_close
    while end + 1 < len(structure) and structure[end + 1] == '.':
        end += 1

    # Slice the sequence and structure
    return {
        "start": start,
        "end": end,
        "hairpin_seq": sequence[start:end + 1],
        "hairpin_struct": structure[start:end + 1]
    }


def trim_trailing_unpaired(sequence, structure):
    """
    Removes nucleotides corresponding to trailing dots in RNA structure.
    """
    # Count trailing dots
    trailing_dots = len(structure) - len(structure.rstrip('.'))

    if trailing_dots == 0:
        return sequence  # nothing to trim
    else:
        return sequence[:-trailing_dots]

def find_rbs_in_hairpin(hairpin_seq):
    """
    Finds the Shine-Dalgarno-like sequence in a terminal hairpin.

    Args:
        hairpin_seq (str): The RNA sequence of the terminal hairpin.

    Returns:
        dict: {
            'found_rbs': bool,       # True if a G-rich 6-mer found 5-13 nt upstream of AUG
            'aug_index': int,        # Index of last AUG in hairpin_seq
            'rbs_seq': str or None,  # The 6-nt G-rich Shine-Dalgarno candidate
            'rbs_region': str        # Full upstream region scanned
        }
    """
    seq = hairpin_seq.upper()
    last_aug = seq.rfind("AUG")
    if last_aug == -1:
        return {
            "found_rbs": False,
            "aug_index": None,
            "rbs_seq": None,
            "rbs_region": None
        }

    # Search 5-13 nt upstream of AUG
    search_start = max(0, last_aug - 13)
    search_end = max(0, last_aug - 5)
    rbs_region = seq[search_start:search_end]

    found = False
    rbs_seq = None
    for i in range(len(rbs_region) - 5):
        window = rbs_region[i:i + 6]
        if window.count("G") >= 3:
            found = True
            rbs_seq = window
            break  # Take the first valid G-rich window

    return {
        "found_rbs": found,
        "aug_index": last_aug,
        "rbs_seq": rbs_seq,
        "rbs_region": rbs_region
    }


def get_rbs_dot_struct(rbs_seq, hairpin_seq, hairpin_struct):
    """
        Extract the RBS dot structure from the hairpin dot structure.

        Args:
            rbs_seq (str): RBS sequence (e.g., "AGGAGG")
            hairpin_seq (str): Full hairpin sequence
            hairpin_struct (str): Dot-bracket structure of the hairpin

        Returns:
            str or None: Dot-bracket structure corresponding to the RBS sequence,
                         or None if RBS not found or structure length mismatch
        """
    # Handle None or empty inputs
    if not rbs_seq or not hairpin_seq or not hairpin_struct:
        return None

    # Validate that structure and sequence lengths match
    if len(hairpin_seq) != len(hairpin_struct):
        return None

    # Find the RBS sequence in the hairpin
    rbs_start = hairpin_seq.find(rbs_seq)

    if rbs_start == -1:
        return None

    # Extract the corresponding structure
    rbs_end = rbs_start + len(rbs_seq)
    rbs_struct = hairpin_struct[rbs_start:rbs_end]

    return rbs_struct

def calc_rbs_paired_percent(rbs_struct):
    """
    Calculate the percentage of paired nucleotides in the RBS structure.

    Args:
        rbs_struct (str): Dot-bracket structure of the RBS sequence
                         ('.' = unpaired, '(' or ')' = paired)

    Returns:
        float: Percentage of paired nucleotides (0-100)
    """
    if len(rbs_struct) == 0:
        return 0.0

    # Count paired positions (both '(' and ')')
    paired_count = sum(1 for char in rbs_struct if char in '()')

    # Calculate percentage
    percent_paired = (paired_count / len(rbs_struct)) * 100

    return percent_paired


def validate_sequence(sequence: str, allowed_chars: str = "ACGU") -> bool:
    """
    Check if sequence contains only valid nucleotides

    Args:
        sequence: RNA/DNA sequence string
        allowed_chars: String of allowed characters (default: "ACGU")

    Returns:
        True if valid, False otherwise
    """
    return all(c in allowed_chars for c in sequence.upper())


def base_pair_percentages(sequence, structure):
    """
    Calculate AU%, GC%, GU% for all paired nucleotides.
    sequence: RNA sequence string
    structure: dot-bracket structure string
    """
    stack = []
    pairs = []

    # Find all paired positions
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            j = stack.pop()
            pairs.append((j, i))

    # Count base pair types
    counts = {'AU': 0, 'UA': 0, 'GC': 0, 'CG': 0, 'GU': 0, 'UG': 0}
    for i, j in pairs:
        pair = sequence[i] + sequence[j]
        if pair in counts:
            counts[pair] += 1

    total_pairs = len(pairs)
    AU_percent = (counts['AU'] + counts['UA']) / total_pairs * 100 if total_pairs else 0
    GC_percent = (counts['GC'] + counts['CG']) / total_pairs * 100 if total_pairs else 0
    GU_percent = (counts['GU'] + counts['UG']) / total_pairs * 100 if total_pairs else 0

    return AU_percent, GC_percent, GU_percent


def gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)

def g_content (seq):
    g = seq.count("G")
    return(g)

def c_content (seq):
    c = seq.count("C")
    return(c)


def mfe_in_range(mfe, min_val, max_val):
    """
    Check if an MFE value is between min_val and max_val.

    Args:
        mfe (float or str): Minimum Free Energy value from RNAfold.
        min_val (float): Lower bound (default -15).
        max_val (float): Upper bound (default -5).

    Returns:
        bool: True if mfe is within range, False otherwise.
    """
    try:
        # Convert to float in case mfe is a string (with or without parentheses)
        mfe_float = float(str(mfe).strip("()"))
        return min_val <= mfe_float <= max_val
    except (ValueError, TypeError):
        # In case conversion fails
        return False

def base_pair_in_range(base_content, min_val, max_val):
    try:
        # Convert to float (strip parentheses if present)
        base_content = float(str(base_content).strip("()"))
        return min_val <= base_content <= max_val
    except (ValueError, TypeError):
        # Return False if conversion fails
        return False




def hairpin_mfe_at_temps(hairpin_seq, temps=[25, 37, 42]):
    mfe_results = {}

    for temp in temps:
        md = RNA.md()  # Create a model details object
        md.temperature = float(temp)
        md.noLP = 1  # ✅ Match RNAfold web server

        fc = RNA.fold_compound(hairpin_seq, md)  # Pass md at creation
        structure, mfe = fc.mfe()
        mfe_results[temp] = (structure, mfe)

    return mfe_results

#
def fold_at_temp(seq, temp):
    md = RNA.md()
    md.temperature = float(temp)
    md.dangles = 2  # Default
    md.noLP = 1  # 0 Allow lonely pairs (default), 1 dont allow LP
    md.noGU = 0  # Allow GU pairs (default)

    fc = RNA.fold_compound(seq, md)
    structure, mfe = fc.mfe()
    return structure, mfe


def base_pairs_at_temps_struct(hairpin_seq, temp=25):
        md = RNA.md()  # Create a model details object
        md.temperature = float(temp)
        md.noLP = 1  # ✅ Match RNAfold web server

        fc = RNA.fold_compound(hairpin_seq, md)  # Pass md at creation
        structure, mfe = fc.mfe()
        base_pair_temp_struct = structure
        return base_pair_temp_struct


def save_results_to_csv(results: List[Dict[str, Any]], output_file: Path,
                        temps: List[int] = [25, 37, 42]) -> None:
    """Save results to CSV file"""
    import csv

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)

        # Dynamic headers based on temperatures
        headers = [
            "Name", "Sequence", "Structure", "MFE",
            "Hairpin Position", "Hairpin Sequence", "Hairpin Structure"
        ]
        headers.extend([f"MFE at {t}°C" for t in temps])
        headers.extend(["PNG Original", "PNG Hairpin"])

        writer.writerow(headers)

        for result in results:
            row = [
                result['name'],
                result['sequence'],
                result['structure'],
                f"{result['mfe']:.2f}",
                f"{result['hairpin_start']}-{result['hairpin_end']}",
                result['hairpin_seq'],
                result['hairpin_struct']
            ]
            row.extend([result.get(f'mfe_{t}', 'N/A') for t in temps])
            row.extend([result['png_original'], result['png_hairpin']])

            writer.writerow(row)

    print(f"✅ Results saved to {output_file}")


def calculate_results_final(
        sequences: List[Tuple[str, str]],
        output_dir: Path,
        settings: Dict[str, int],
        progress_callback: Optional[Callable[[str], None]] = None,
        csv_settings_manager = None  # ✨ NEW: Accept settings manager
) -> List[Dict[str, Any]]:
    """
    Analyze RNA sequences for thermometer properties

    Args:
        sequences: List of (name, sequence) tuples
        output_dir: Directory for output files
        progress_callback: Optional function to call with progress messages

    Returns:
        List of result tuples
    """
    results = []

    def log(message: str):
        """Helper to log messages to both console and GUI"""
        print(message)
        if progress_callback:
            progress_callback(message)

    log(f"🧬 Analyzing {len(sequences)} RNA sequences...\n")

    # Create structures subdirectory
    structures_dir = output_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    total = len(sequences)
    #total = 15

    for idx, (og_name, og_seq) in enumerate(sequences, 1):
        log(f"\n{'=' * 60}")
        log(f"[{idx}/{total}] Processing: {og_name}")
        log(f"  Sequence length: {len(og_seq)} nt")
        found_count = 0;

        # Skip very short sequences
        if len(og_seq) <= 4:
            log(f"  Sequence too short for hairpin detection, skipping.\n")
            continue

            # ✨ NEW: Apply sequence preprocessing
        if csv_settings_manager:
            seq_settings = csv_settings_manager.settings.get("sequence_processing", {})
            if seq_settings.get("append_sequence_enabled", False):
                append_seq = seq_settings.get("append_sequence", "AUG").upper()
                position = seq_settings.get("append_position", "end")

                if position == "start":
                    og_seq = append_seq + og_seq
                    log(f"  ✨ Prepended '{append_seq}' to sequence (5' end)")
                else:
                    og_seq = og_seq + append_seq
                    log(f"  ✨ Appended '{append_seq}' to sequence (3' end)")

                log(f"  Modified sequence length: {len(og_seq)} nt")





        # Get calculation settings
        calc_settings = {}
        if csv_settings_manager:
            calc_settings = csv_settings_manager.settings.get("calculation_settings", {})



        # ✨ NEW: Calculate composition for ORIGINAL sequence
        original_comp = {"AU%": 0, "GC%": 0, "GU%": 0}
        if calc_settings.get("calculate_original_composition", False):
            log(f"  Calculating original sequence composition...")
            original_comp = calculate_composition(og_seq)

            # ✨ CONDITIONAL: Original sequence MFE at temps (only if needed)
        mfe_25_og = mfe_37_og = mfe_42_og = 0.0
        structure_25 = structure_37 = structure_42 = ""
        if calc_settings.get("calculate_original_mfe_temps", False):
            log(f"  Folding original sequence at 25°C, 37°C, 42°C...")
            structure_25, mfe_25_og = fold_at_temp(og_seq, 25)
            structure_37, mfe_37_og = fold_at_temp(og_seq, 37)
            structure_42, mfe_42_og = fold_at_temp(og_seq, 42)
        else:
            # Still need structure at 25°C for hairpin detection
            log(f"  Folding at 25°C (for hairpin detection)...")
            structure_25, mfe_25_og = fold_at_temp(og_seq, 25)

        log(f"  MFE: {mfe_25_og:.2f} kcal/mol")


        # ✨ NEW: Check if original sequence values are in range
        orig_mfe_25_in_range = mfe_in_range(mfe_25_og, settings.get('orig_mfe_25_min', -100),
                                            settings.get('orig_mfe_25_max', 100))
        orig_mfe_37_in_range = mfe_in_range(mfe_37_og, settings.get('orig_mfe_37_min', -100),
                                            settings.get('orig_mfe_37_max', 100))
        orig_mfe_42_in_range = mfe_in_range(mfe_42_og, settings.get('orig_mfe_42_min', -100),
                                            settings.get('orig_mfe_42_max', 100))

        orig_au_in_range = base_pair_in_range(original_comp["AU%"], settings.get('orig_au_min', 0),
                                              settings.get('orig_au_max', 100))
        orig_gc_in_range = base_pair_in_range(original_comp["GC%"], settings.get('orig_gc_min', 0),
                                              settings.get('orig_gc_max', 100))
        orig_gu_in_range = base_pair_in_range(original_comp["GU%"], settings.get('orig_gu_min', 0),
                                              settings.get('orig_gu_max', 100))

        # Convert to strings
        orig_mfe_25_str = "In Range" if orig_mfe_25_in_range else "Not in Range"
        orig_mfe_37_str = "In Range" if orig_mfe_37_in_range else "Not in Range"
        orig_mfe_42_str = "In Range" if orig_mfe_42_in_range else "Not in Range"
        orig_au_str = "In Range" if orig_au_in_range else "Not in Range"
        orig_gc_str = "In Range" if orig_gc_in_range else "Not in Range"
        orig_gu_str = "In Range" if orig_gu_in_range else "Not in Range"

        # ✨ NEW: Calculate original sequence quality score (0-6)
        orig_quality_score = sum([
            orig_mfe_25_in_range,
            orig_mfe_37_in_range,
            orig_mfe_42_in_range,
            orig_au_in_range,
            orig_gc_in_range,
            orig_gu_in_range
        ])

        # Log original sequence filter results
        log(f"  Original sequence filters:")
        log(f"    MFE 25°C: {'✓' if orig_mfe_25_in_range else '✗'} {orig_mfe_25_str}")
        log(f"    MFE 37°C: {'✓' if orig_mfe_37_in_range else '✗'} {orig_mfe_37_str}")
        log(f"    MFE 42°C: {'✓' if orig_mfe_42_in_range else '✗'} {orig_mfe_42_str}")
        log(f"    AU%: {'✓' if orig_au_in_range else '✗'} {orig_au_str}")
        log(f"    GC%: {'✓' if orig_gc_in_range else '✗'} {orig_gc_str}")
        log(f"    GU%: {'✓' if orig_gu_in_range else '✗'} {orig_gu_str}")
        log(f"  Original Quality Score: {orig_quality_score}/6")


        # Terminal Hairpin Info
        log(f"  Detecting terminal hairpin...")
        term_results = get_terminal_hairpin_with_tail(og_seq, structure_25)


        # Check if a hairpin was detected
        if term_results is None or term_results.get("hairpin_seq") is None:
            log(f"  No terminal hairpin detected, skipping this sequence.\n")
            continue  # skip to next sequence
        hairpin_seq = term_results["hairpin_seq"]
        hairpin_struct = term_results["hairpin_struct"] # This is the structure from 25°C fold
        hairpin_seq_trimmed = trim_trailing_unpaired(hairpin_seq, hairpin_struct)

       #log(f"  Terminal hairpin: position {start}-{end}")
        log(f"  Hairpin length: {len(hairpin_seq)} nt (trimmed: {len(hairpin_seq_trimmed)} nt)")

        # RBS region
        # ✨ CONDITIONAL: RBS region (only if enabled)
        RBS_seq = None
        RBS_dot_struct = None
        RBS_paired_percent = None

        if calc_settings.get("calculate_rbs", True):
            log(f"  Searching for RBS...")
            RBS_results = find_rbs_in_hairpin(hairpin_seq)
            RBS_seq = RBS_results["rbs_seq"]

            if RBS_seq:
                log(f"  ✓ RBS found: {RBS_seq}")
                RBS_dot_struct = get_rbs_dot_struct(RBS_seq, hairpin_seq, hairpin_struct)
                if RBS_dot_struct is not None:
                    RBS_paired_percent = calc_rbs_paired_percent(RBS_dot_struct)
                    log(f"  RBS structure: {RBS_dot_struct}")
                    log(f"  RBS paired: {RBS_paired_percent:.1f}%")
            else:
                log(f"  ✗ No RBS detected")
        else:
            log(f"  Skipping RBS calculation (disabled in settings)")



        # MFE at different temperatures
        log(f"  Calculating MFE at 25°C, 37°C, 42°C...")
        MFE_results = hairpin_mfe_at_temps(hairpin_seq_trimmed, temps=[25, 37, 42])

        # Extract MFE values
        mfe_25 = MFE_results[25][1]
        mfe_37 = MFE_results[37][1]
        mfe_42 = MFE_results[42][1]

        # Check if in range
        mfe_25_in_range = mfe_in_range(mfe_25, settings['mfe_25_min'],  settings['mfe_25_max'])
        mfe_37_in_range = mfe_in_range(mfe_37,  settings['mfe_37_min'],  settings['mfe_37_max'])
        mfe_42_in_range = mfe_in_range(mfe_42, settings['mfe_42_min'],settings['mfe_42_max'])

        mfe_25_str = "In Range" if mfe_25_in_range else "Not in Range"
        if mfe_25_str == "In Range":
            found_count+=1
        mfe_37_str = "In Range" if mfe_37_in_range else "Not in Range"
        if mfe_37_str == "In Range":
            found_count += 1
        mfe_42_str = "In Range" if mfe_42_in_range else "Not in Range"
        if mfe_42_str == "In Range":
            found_count += 1

        status_25 = "✓" if mfe_25_in_range else "✗"
        status_37 = "✓" if mfe_37_in_range else "✗"
        status_42 = "✓" if mfe_42_in_range else "✗"

        log(f"    25°C: {mfe_25:6.2f} kcal/mol {status_25} {mfe_25_str}")
        log(f"    37°C: {mfe_37:6.2f} kcal/mol {status_37} {mfe_37_str}")
        log(f"    42°C: {mfe_42:6.2f} kcal/mol {status_42} {mfe_42_str}")


        # Base pair composition - use the ORIGINAL hairpin structure from 25°C
        log(f"  Analyzing base pair composition...")
        AU, GC, GU = base_pair_percentages(hairpin_seq, hairpin_struct)  # ✅ Use original structure

        # Check if in range
        AU_in_range = base_pair_in_range(AU, settings['au_min'], settings['au_max'])
        GC_in_range = base_pair_in_range(GC, settings['gc_min'], settings['gc_max'])
        GU_in_range = base_pair_in_range(GU, settings['gu_min'], settings['gu_max'])

        AU_str = "In Range" if AU_in_range else "Not in Range"
        if AU_str == "In Range":
            found_count += 1
        GC_str = "In Range" if GC_in_range else "Not in Range"
        if GC_str == "In Range":
            found_count += 1
        GU_str = "In Range" if GU_in_range else "Not in Range"
        if GU_str == "In Range":
            found_count += 1

        status_AU = "✓" if AU_in_range else "✗"
        status_GC = "✓" if GC_in_range else "✗"
        status_GU = "✓" if GU_in_range else "✗"

        log(f"    AU: {AU:5.1f}% {status_AU} {AU_str}")
        log(f"    GC: {GC:5.1f}% {status_GC} {GC_str}")
        log(f"    GU: {GU:5.1f}% {status_GU} {GU_str}")

        # Generate structure diagrams
        #log(f"  Generating structure diagrams...")
        #hyperlink_original = ""
        #hyperlink_hairpin = ""

        total_found_count = found_count


        # ✨ CHANGED: Store as dictionary instead of tuple
        # THIS IS WHERE RESULTS ARE SHOWN, ADD HERE IF NEW RESULTS ADDED
        result_data = {
            "name": og_name,
            "original_sequence": og_seq,
            "original_structure": structure_25,
            "original_mfe_25": f"{mfe_25_og:.2f}",
            "original_mfe_37": f"{mfe_37_og:.2f}",
            "original_mfe_42": f"{mfe_42_og:.2f}",
            "original_au_percent": original_comp["AU%"],  # 🔑 NEW
            "original_gc_percent": original_comp["GC%"],  # 🔑 NEW
            "original_gu_percent": original_comp["GU%"],  # 🔑 NEW

            # ✨ NEW: Original sequence range checks
            "original_mfe_25_in_range": orig_mfe_25_str,
            "original_mfe_37_in_range": orig_mfe_37_str,
            "original_mfe_42_in_range": orig_mfe_42_str,
            "original_au_in_range": orig_au_str,
            "original_gc_in_range": orig_gc_str,
            "original_gu_in_range": orig_gu_str,

            "hairpin_sequence": hairpin_seq,
            "hairpin_structure": hairpin_struct,
            "hairpin_au_percent": AU,  # Already calculated
            "hairpin_gc_percent": GC,
            "hairpin_gu_percent": GU,
            "mfe_25c_hairpin": f"{mfe_25:.2f}",
            "mfe_37c_hairpin": f"{mfe_37:.2f}",
            "mfe_42c_hairpin": f"{mfe_42:.2f}",
            "mfe_25_in_range_hairpin": mfe_25_str,
            "mfe_37_in_range_hairpin": mfe_37_str,
            "mfe_42_in_range_hairpin": mfe_42_str,
            "au_in_range_hairpin": AU_str,
            "gc_in_range_hairpin": GC_str,
            "gu_in_range_hairpin": GU_str,
            "rbs_sequence": RBS_seq if RBS_seq else "Not Found",
            "rbs_structure": RBS_dot_struct if RBS_dot_struct else "N/A",
            "rbs_paired_percent": f"{RBS_paired_percent:.2f}" if RBS_paired_percent is not None else "N/A",
            "quality_score_hairpin": total_found_count,
            "quality_score_original": orig_quality_score  # ✨ NEW
        }

        results.append(result_data)

        log(f"  ✓ Completed {og_name} ({idx}/{total})")

     # ===== SORT RESULTS (BEFORE SAVING) =====
    log(f"\n{'=' * 60}")
    log(f"📊 Sorting results by Total In Range Count...")
    results.sort(key=lambda x: x.get("quality_score_hairpin", 0), reverse=True)
    log(f"✅ Sorted! Best candidates at top.")

    # ===== HIGHLIGHT TOP CANDIDATES =====
    top_candidates = [r for r in results if r.get("quality_score_hairpin", 0) >= 4]
    if top_candidates:
        log(f"\n🎯 Top Candidates (4+ criteria): {len(top_candidates)} sequences")
        for result in top_candidates[:5]:  # Show first 5
            log(f"   • {result['name']} - {result.get('quality_score_hairpin', 0)}/6 criteria")
        if len(top_candidates) > 5:
            log(f"   ... and {len(top_candidates) - 5} more")

        # Save results to CSV using settings
    log(f"\n{'=' * 60}")
    log(f"💾 Saving results to CSV...")

    output_file = output_dir / "rna_results.csv"

    # ✨ NEW: Load CSV settings
    try:
        csv_settings = SettingsManager("csv_output_settings.json")
        headers = csv_settings.get_enabled_columns()
        log(f"📊 Using custom CSV columns: {len(headers)} columns")
    except:
        # Fallback to all columns if settings not available
        log("⚠ Using default CSV columns (all columns)")
        headers = [
            "Name",
            "Sequence",
            "Structure",
            "Original_MFE_25C",
            "Original_MFE_37C",
            "Original_MFE_42C",
            "Original_AU%",
            "Original_GC%",
            "Original_GU%",
            "Original_MFE_25C_InRange",
            "Original_MFE_37C_InRange",
            "Original_MFE_42C_InRange",
            "Original_AU%_InRange",
            "Original_GC%_InRange",
            "Original_GU%_InRange",
            "Hairpin_Sequence",
            "Hairpin_Structure",
            "Hairpin_AU%",
            "Hairpin_GC%",
            "Hairpin_GU%",
            "Hairpin_MFE_25C",
            "Hairpin_MFE_37C",
            "Hairpin_MFE_42C",
            "Hairpin_MFE_25C_InRange",
            "Hairpin_MFE_37C_InRange",
            "Hairpin_MFE_42C_InRange",
            "Hairpin_AU%_InRange",
            "Hairpin_GC%_InRange",
            "Hairpin_GU%_InRange",
            "RBS_Sequence",
            "RBS_Structure",
            "RBS_Paired%",
            "Quality_Score_Hairpin",
            "Quality_Score_Original"
        ]
        csv_settings = None

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        # Write data rows
        for result_data in results:
            if csv_settings:
                # Use settings to build row
                from RnaThermofinder.utils.analysis_helpers import build_csv_row
                row = build_csv_row(result_data, csv_settings)

                # Debug: Check if row is valid
                if row is None or not isinstance(row, list):
                    log(f"⚠ Warning: build_csv_row returned invalid data, using fallback")
                    csv_settings = None  # Switch to fallback mode

            if not csv_settings:
                # Fallback: write all columns
                row = [
                    result_data.get("name", ""),
                    result_data.get("original_sequence", ""),
                    result_data.get("original_structure", ""),
                    result_data.get("original_mfe_25", ""),
                    result_data.get("original_mfe_37", ""),
                    result_data.get("original_mfe_42", ""),
                    result_data.get("original_au_percent", ""),
                    result_data.get("original_gc_percent", ""),
                    result_data.get("original_gu_percent", ""),

                    result_data.get("original_mfe_25_in_range", ""),
                    result_data.get("original_mfe_37_in_range", ""),
                    result_data.get("original_mfe_42_in_range", ""),
                    result_data.get("original_au_in_range", ""),
                    result_data.get("original_gc_in_range", ""),
                    result_data.get("original_gu_in_range", ""),

                    result_data.get("hairpin_sequence", ""),
                    result_data.get("hairpin_structure", ""),
                    result_data.get("hairpin_au_percent", ""),
                    result_data.get("hairpin_gc_percent", ""),
                    result_data.get("hairpin_gu_percent", ""),
                    result_data.get("mfe_25c_hairpin", ""),
                    result_data.get("mfe_37c_hairpin", ""),
                    result_data.get("mfe_42c_hairpin", ""),
                    result_data.get("mfe_25_in_range_hairpin", ""),
                    result_data.get("mfe_37_in_range_hairpin", ""),
                    result_data.get("mfe_42_in_range_hairpin", ""),
                    result_data.get("au_in_range_hairpin", ""),
                    result_data.get("gc_in_range_hairpin", ""),
                    result_data.get("gu_in_range_hairpin", ""),
                    result_data.get("rbs_sequence", ""),
                    result_data.get("rbs_structure", ""),
                    result_data.get("rbs_paired_percent", ""),
                    result_data.get("quality_score_hairpin", ""),
                    result_data.get("quality_score_original", ""),
                ]

            writer.writerow(row)
    log(f"✅ All results saved to: {output_file.name}")
    log(f"✅ Analysis complete! Processed {len(results)} sequences")

    return results
from pathlib import Path

import RNA
from typing import List, Tuple, Callable, Optional, Dict, Any
import csv

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

        fc = RNA.fold_compound(hairpin_seq, md)  # Pass md at creation
        structure, mfe = fc.mfe()
        mfe_results[temp] = (structure, mfe)

    return mfe_results

def fold_at_temps(seq, temp):
    md = RNA.md()  # Create a model details object
    md.temperature = float(temp)
    fc = RNA.fold_compound(seq, md)  # Pass md at creation
    structure, mfe = fc.mfe()
    return structure,mfe



def base_pairs_at_temps_struct(hairpin_seq, temp=25):
        md = RNA.md()  # Create a model details object
        md.temperature = float(temp)

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

    for idx, (og_name, og_seq) in enumerate(sequences, 1):
        log(f"\n{'=' * 60}")
        log(f"[{idx}/{total}] Processing: {og_name}")
        log(f"  Sequence length: {len(og_seq)} nt")
        found_count = 0;

        # Skip very short sequences
        if len(og_seq) <= 4:
            log(f"  Sequence too short for hairpin detection, skipping.\n")
            continue

        # Add AUG at the end
        og_seq = og_seq + "AUG"


        # Sequence fold info at 25°C
        log(f"  Folding at 25°C...")
        structure, mfe = fold_at_temps(og_seq, 25)
        log(f"  MFE: {mfe:.2f} kcal/mol")

        # Terminal Hairpin Info
        log(f"  Detecting terminal hairpin...")
        term_results = get_terminal_hairpin_with_tail(og_seq, structure)
        # Check if a hairpin was detected
        if term_results is None or term_results.get("hairpin_seq") is None:
            log(f"  No terminal hairpin detected, skipping this sequence.\n")
            continue  # skip to next sequence
        #start = term_results["start"]
        #end = term_results["end"]
        hairpin_seq = term_results["hairpin_seq"]
        hairpin_struct = term_results["hairpin_struct"]
        hairpin_seq_trimmed = trim_trailing_unpaired(hairpin_seq, hairpin_struct)

       # log(f"  Terminal hairpin: position {start}-{end}")
        log(f"  Hairpin length: {len(hairpin_seq)} nt (trimmed: {len(hairpin_seq_trimmed)} nt)")

        # RBS region
        log(f"  Searching for RBS...")
        RBS_results = find_rbs_in_hairpin(hairpin_seq)
        RBS_seq = RBS_results["rbs_seq"]

        # Get RBS dot structure and calculate pairing percentage
        RBS_dot_struct = None
        RBS_paired_percent = None

        if RBS_seq:
            log(f"  ✓ RBS found: {RBS_seq}")

            # Extract RBS structure from hairpin
            RBS_dot_struct = get_rbs_dot_struct(RBS_seq, hairpin_seq, hairpin_struct)

            if RBS_dot_struct is not None:
                # Calculate percentage of paired nucleotides
                RBS_paired_percent = calc_rbs_paired_percent(RBS_dot_struct)
                log(f"  RBS structure: {RBS_dot_struct}")
                log(f"  RBS paired: {RBS_paired_percent:.1f}%")
            else:
                log(f"  ⚠ Could not extract RBS structure")
        else:
            log(f"  ✗ No RBS detected")

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

        # Base pair composition
        log(f"  Analyzing base pair composition...")
        base_pair_temp_struct_hairpin = base_pairs_at_temps_struct(hairpin_seq, 25)
        AU, GC, GU = base_pair_percentages(hairpin_seq, base_pair_temp_struct_hairpin)

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
        log(f"  Generating structure diagrams...")

        # Save original structure SVG
        #original_svg = structures_dir / f"{og_name}_structure.svg"
        #RNA.svg_rna_plot(og_seq, structure, str(original_svg))

        # Save hairpin structure SVG (fixed filename - was same as original)
        #hairpin_svg = structures_dir / f"{og_name}_terminalhairpin_structure.svg"
       # RNA.svg_rna_plot(hairpin_seq, hairpin_struct, str(hairpin_svg))

       # log(f"  ✓ Saved: {original_svg.name}")
       # log(f"  ✓ Saved: {hairpin_svg.name}")

        # Create hyperlinks for Excel (with proper encoding)
        import urllib.parse

       # svg_og_path = original_svg.resolve()
       # encoded_og_path = urllib.parse.quote(svg_og_path.as_posix())
       # file_og_uri = f"file://{encoded_og_path}"
        #hyperlink_original = f'=HYPERLINK("{file_og_uri}", "View Structure")'
        hyperlink_original = ""

        #svg_pin_path = hairpin_svg.resolve()
       # encoded_pin_path = urllib.parse.quote(svg_pin_path.as_posix())
       # file_pin_uri = f"file://{encoded_pin_path}"
       # hyperlink_hairpin = f'=HYPERLINK("{file_pin_uri}", "View Structure")'
        hyperlink_hairpin = ""

        total_found_count = found_count

        # Append all data as a single tuple
        results.append((
            og_name,
            og_seq,
            structure,
            hairpin_seq,
            hairpin_struct,
            RBS_seq if RBS_seq else "Not Found",
            RBS_dot_struct if RBS_dot_struct else "N/A",
            f"{RBS_paired_percent:.2f}" if RBS_paired_percent is not None else "N/A",
            f"{mfe_25:.2f}",
            f"{mfe_37:.2f}",
            f"{mfe_42:.2f}",
            mfe_25_str,
            mfe_37_str,
            mfe_42_str,
            f"{AU:.2f}",
            f"{GC:.2f}",
            f"{GU:.2f}",
            AU_str,
            GC_str,
            GU_str,
            hyperlink_original,
            hyperlink_hairpin,
            total_found_count

        ))

        log(f"  ✓ Completed {og_name} ({idx}/{total})")

     # ===== SORT RESULTS (BEFORE SAVING) =====
    log(f"\n{'=' * 60}")
    log(f"📊 Sorting results by Total In Range Count...")
    results.sort(key=lambda x: x[-1], reverse=True)
    log(f"✅ Sorted! Best candidates at top.")

    # ===== HIGHLIGHT TOP CANDIDATES =====
    top_candidates = [r for r in results if r[-1] >= 4]
    if top_candidates:
        log(f"\n🎯 Top Candidates (4+ criteria): {len(top_candidates)} sequences")
        for result in top_candidates[:5]:  # Show first 5
            log(f"   • {result[0]} - {result[-1]}/6 criteria")
        if len(top_candidates) > 5:
            log(f"   ... and {len(top_candidates) - 5} more")



    # Save results to CSV
    log(f"\n{'=' * 60}")
    log(f"💾 Saving results to CSV...")

    output_file = output_dir / "rna_results.csv"
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow([
            "Name",
            "Sequence",
            "Structure",
            "Hairpin Sequence",
            "Hairpin Structure",
            "RBS Sequence (6nt window)",
            "RBS Dot Structure",
            "RBS Paired %",
            "MFE at 25C (kcal/mol)",
            "MFE at 37C (kcal/mol)",
            "MFE at 42C (kcal/mol)",
            "MFE in Range at 25C",
            "MFE in Range at 37C",
            "MFE in Range at 42C",
            "AU% Composition",
            "GC% Composition",
            "GU% Composition",
            "AU% in Range",
            "GC% in Range",
            "GU% in Range",
            "SVG Original",
            "SVG Hairpin",
            "Total In Range Count"
        ])

        # Write rows
        for row in results:
            writer.writerow(row)

    log(f"✅ All results saved to: {output_file.name}")
    log(f"✅ Analysis complete! Processed {len(results)} sequences")

    return results
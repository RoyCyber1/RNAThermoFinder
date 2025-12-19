from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
import os
import time

Entrez.email = "roycyber14@gmail.com"

# Your sequence OR try E. coli sequence for testing
my_sequence = "AACAGTGTAAAAACATTGTTTTTCCGGTCTTTTTTTGTGTGCCGGAAGGAGAGTGGAAAGGATG"
result_file = "blast_results.xml"

if os.path.exists(result_file):
    print("Using saved BLAST results...")
    result_handle = open(result_file)
else:
    print("Running BLAST... (this may take 2-5 minutes)")
    start_time = time.time()

    result_handle = NCBIWWW.qblast(
        "blastn",
        "nt",  # Databases "nt", "refseq_genomic"
        my_sequence,
        hitlist_size=10,  # More hits to find annotated ones
        expect=0.01
    )

    elapsed = time.time() - start_time
    print(f"✓ BLAST completed in {elapsed:.1f} seconds")

    blast_output = result_handle.read()
    with open(result_file, "w") as out:
        out.write(blast_output)
    print("Results saved to", result_file)

    result_handle = open(result_file)

print("\n" + "=" * 60)
print("Parsing BLAST results...")
print("=" * 60)

blast_records = NCBIXML.parse(result_handle)

hit_count = 0
annotated_hits = 0

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.01:
                hit_count += 1
                accession = alignment.accession

                print(f"\n{'=' * 60}")
                print(f"HIT #{hit_count}")
                print(f"{'=' * 60}")
                print(f"Match: {alignment.title}")
                print(f"Accession: {accession}")
                print(f"E-value: {hsp.expect}")
                print(f"Identity: {hsp.identities}/{hsp.align_length} ({100 * hsp.identities / hsp.align_length:.1f}%)")

                print("\nFetching GenBank record...")
                try:
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=accession,
                        rettype="gb",
                        retmode="text"
                    )

                    record = SeqIO.read(fetch_handle, "genbank")
                    print(f"✓ Record ID: {record.id}")
                    print(f"✓ Description: {record.description}")
                    print(f"✓ Sequence length: {len(record.seq):,} bp")

                    hit_start = min(hsp.sbjct_start, hsp.sbjct_end)
                    hit_end = max(hsp.sbjct_start, hsp.sbjct_end)

                    print(f"\nMatch position: {hit_start:,} - {hit_end:,}")
                    print("\n" + "-" * 60)
                    print("GENOMIC CONTEXT ANALYSIS:")
                    print("-" * 60)

                    # Try to find gene OR CDS features
                    features_list = []
                    for feature in record.features:
                        if feature.type in ["gene", "CDS"]:
                            start = int(feature.location.start)
                            end = int(feature.location.end)

                            # Try different qualifiers for name
                            name = (feature.qualifiers.get('gene',
                                                           feature.qualifiers.get('product',
                                                                                  feature.qualifiers.get('locus_tag',
                                                                                                         ['Unknown'])))[
                                0])

                            features_list.append((start, end, name, feature.type))

                    if not features_list:
                        print("⚠ No gene or CDS annotations found in this record")
                        print("   (This is common with newly sequenced genomes)")
                    else:
                        annotated_hits += 1
                        features_list.sort()

                        position_found = False

                        for feat_start, feat_end, feat_name, feat_type in features_list:
                            # INSIDE
                            if hit_start >= feat_start and hit_end <= feat_end:
                                print(f"\n🎯 INSIDE {feat_type.upper()}: '{feat_name}'")
                                print(f"   {feat_type} position: {feat_start:,} - {feat_end:,}")
                                print(f"   Match position: {hit_start:,} - {hit_end:,}")
                                print(f"   Distance from {feat_type} start: {hit_start - feat_start:,} bp")
                                print(f"   Distance from {feat_type} end: {feat_end - hit_end:,} bp")
                                position_found = True
                                break

                            # BEFORE - match ends before feature starts
                            elif hit_end < feat_start:
                                distance = feat_start - hit_end
                                print(f"\n⬅️  BEFORE {feat_type.upper()}: '{feat_name}'")
                                print(f"   {feat_type} position: {feat_start:,} - {feat_end:,}")
                                print(f"   Match position: {hit_start:,} - {hit_end:,}")
                                print(f"   Distance to {feat_type}: {distance:,} bp upstream")
                                position_found = True
                                break

                        # AFTER - Check if after the last feature OR between features
                        if not position_found:
                            # Check if we're after ANY feature
                            for feat_start, feat_end, feat_name, feat_type in features_list:
                                # Match starts after this feature ends
                                if hit_start > feat_end:
                                    # Check if there's a next feature
                                    next_feature_exists = False
                                    for next_start, next_end, next_name, next_type in features_list:
                                        if next_start > feat_end and hit_end < next_start:
                                            # We're in an intergenic region
                                            dist_from = hit_start - feat_end
                                            dist_to = next_start - hit_end
                                            print(f"\n🔀 INTERGENIC REGION")
                                            print(f"   After {feat_type}: '{feat_name}' ({dist_from:,} bp downstream)")
                                            print(f"   Before {next_type}: '{next_name}' ({dist_to:,} bp upstream)")
                                            print(f"   Match position: {hit_start:,} - {hit_end:,}")
                                            position_found = True
                                            next_feature_exists = True
                                            break

                                    if position_found:
                                        break

                                    # If no next feature found and this is the last one we checked
                                    if not next_feature_exists and feat_start == features_list[-1][0]:
                                        distance = hit_start - feat_end
                                        print(f"\n➡️  AFTER {feat_type.upper()}: '{feat_name}'")
                                        print(f"   {feat_type} position: {feat_start:,} - {feat_end:,}")
                                        print(f"   Match position: {hit_start:,} - {hit_end:,}")
                                        print(f"   Distance from {feat_type}: {distance:,} bp downstream")
                                        position_found = True
                                        break

                    print("-" * 60)
                    fetch_handle.close()

                except Exception as e:
                    print(f"✗ Error: {e}")

                time.sleep(0.5)

result_handle.close()

print(f"\n{'=' * 60}")
print(f"Analysis complete!")
print(f"Total hits found: {hit_count}")
print(f"Hits with gene annotations: {annotated_hits}")
if annotated_hits == 0:
    print("\n💡 TIP: Try using a sequence from a well-studied organism")
    print("   (like E. coli, human, mouse) for better annotations")
print(f"{'=' * 60}")
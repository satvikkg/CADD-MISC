#!/usr/bin/env python
"""
Schrödinger MAE/SDF Preprocessor Workflow with Stereochemistry & Conformer Tagging

DESCRIPTION:
Processes MAE or SDF files from Schrödinger Maestro's LigPrep and produces deterministic,
idempotent compound titles that include stereochemistry and conformer indices.

This version uses Schrödinger's native stereochemistry detection and can handle both
MAE and SDF files as input and output.

USAGE:
$SCHRODINGER/run python mae_preprocessor_workflow.py input.mae -o output_basename
$SCHRODINGER/run python mae_preprocessor_workflow.py input.sdf -o output.mae --mae-output

DEPENDENCIES:
- Schrödinger Python API
- pandas (optional, for CSV export)
"""

import argparse
import os
import sys
import gzip
import shutil
import re
import tempfile
from collections import defaultdict

# Schrödinger imports
try:
    from schrodinger.structure import StructureReader, StructureWriter
    from schrodinger.comparison.chirality import get_chirality, get_chiral_centers
    from schrodinger.structutils import analyze
    from schrodinger.utils import log
    SCHRODINGER_AVAILABLE = True
except ImportError:
    SCHRODINGER_AVAILABLE = False
    print("Error: Schrödinger modules not available. Please run with $SCHRODINGER/run python")
    sys.exit(1)

# Optional imports
try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    print("Warning: pandas not available, CSV export disabled")


class MAEPreprocessor:
    """Main class for MAE/SDF preprocessing workflow using Schrödinger API"""
    
    def __init__(self, input_file, output_basename):
        self.input_file = input_file
        self.output_basename = output_basename
        
        # Regex for detecting already tagged compounds
        self.already_tagged_re = re.compile(r'(?:_\d+[RS])+(?:_CONF\d+)?$', re.IGNORECASE)
        
    def read_structures(self):
        """Read structures from MAE or SDF file using Schrödinger StructureReader"""
        structures = []
        
        try:
            with StructureReader(self.input_file) as reader:
                for st in reader:
                    if st is not None:
                        structures.append(st)
            
            print(f"📖 Read {len(structures)} structures from {self.input_file}")
            return structures
            
        except Exception as e:
            print(f"❌ Error reading file {self.input_file}: {str(e)}")
            return []
    
    def normalize_titles(self, structures):
        """Normalize titles: trim whitespace and convert to uppercase"""
        for st in structures:
            if st.title:
                st.title = st.title.strip().upper()
            else:
                st.title = "UNKNOWN"
    
    def group_by_title(self, structures):
        """Group structures by their current title"""
        grouped = defaultdict(list)
        
        for st in structures:
            title = st.title if st.title else "UNKNOWN"
            grouped[title].append(st)
        
        return grouped
    
    def compute_stereochemistry_schrodinger(self, st):
        """Compute stereochemistry signature using Schrödinger's native tools
        
        This function can handle molecules with MULTIPLE chiral centers.
        It returns ALL R/S assignments for all chiral centers in the molecule.
        """
        try:
            # Get comprehensive stereochemistry labels for all atoms
            chirality_labels = get_chirality(st, ignore_annotations=True)
            
            # Get indices of all chiral centers in the molecule
            chiral_centers = get_chiral_centers(st)
            
            print(f"    🔬 Found {len(chiral_centers) if chiral_centers else 0} chiral centers: {chiral_centers}")
            
            # Build stereochemistry signature from ALL chiral centers
            centers = []
            
            if chiral_centers and chirality_labels:
                for atom_idx in chiral_centers:
                    # Convert to 0-based indexing for array access (Schrödinger uses 1-based)
                    zero_based_idx = atom_idx - 1
                    
                    if zero_based_idx < len(chirality_labels):
                        label = chirality_labels[zero_based_idx]
                        
                        # Only include R/S labels (ignore E/Z double bond labels)
                        if label and label.upper() in ['R', 'S']:
                            centers.append((atom_idx, label.upper()))  # Keep 1-based for readability
                            print(f"      📍 Atom {atom_idx}: {label.upper()}")
            
            # Sort by atom index for deterministic, reproducible ordering
            # This ensures "_4R_8S_12R" always comes out the same way
            sorted_centers = sorted(centers, key=lambda x: x[0])
            
            if sorted_centers:
                stereo_string = "_".join([f"{idx}{label}" for idx, label in sorted_centers])
                print(f"    🏷️ Complete stereo signature: {stereo_string}")
            
            return sorted_centers
        
        except Exception as e:
            print(f"Warning: Schrödinger stereochemistry detection failed for structure: {e}")
            return []
    
    def generate_new_titles(self, grouped_structures):
        """Generate new titles with stereochemistry and conformer tags
        
        This handles molecules with MULTIPLE chiral centers by creating
        comprehensive stereochemistry signatures like "_4R_8S_12R_15S"
        """
        new_title_map = {}
        
        for title, group in grouped_structures.items():
            # Skip if already tagged
            if self.already_tagged_re.search(title):
                print(f"⏭️ Skipping already tagged group: '{title}'")
                for st in group:
                    new_title_map[id(st)] = title
                continue
            
            print(f"🔍 Processing group '{title}' with {len(group)} structures")
            
            # Compute stereochemistry for each structure
            stereo_tags = []
            for i, st in enumerate(group):
                print(f"  📋 Analyzing structure {i+1}/{len(group)}")
                centers = self.compute_stereochemistry_schrodinger(st)
                stereo_tags.append(tuple(centers))
            
            # Group by stereo signature (structures with identical stereochemistry)
            stereo_to_structures = defaultdict(list)
            for i, tag in enumerate(stereo_tags):
                stereo_to_structures[tag].append(group[i])
            
            print(f"  🧬 Found {len(stereo_to_structures)} unique stereoisomers in this group")
            
            # Generate final titles
            for stereo_tag, structure_list in stereo_to_structures.items():
                has_stereo = len(stereo_tag) > 0
                
                if has_stereo:
                    # Create comprehensive suffix for ALL chiral centers
                    # Example: [(4, 'R'), (8, 'S'), (12, 'R')] -> "_4R_8S_12R"
                    suffix = "_".join([f"{idx}{conf}" for idx, conf in stereo_tag])
                    print(f"  🧪 Stereo signature with {len(stereo_tag)} centers: {suffix}")
                else:
                    print(f"  ⚪ No stereochemistry detected for this group")
                
                for i, st in enumerate(structure_list, start=1):
                    if has_stereo:
                        if len(structure_list) > 1:
                            # Multiple conformers of the same stereoisomer
                            new_title = f"{title}_{suffix}_CONF{i}"
                        else:
                            # Single structure with this stereochemistry
                            new_title = f"{title}_{suffix}"
                    else:
                        if len(structure_list) > 1:
                            # Multiple conformers, no stereochemistry
                            new_title = f"{title}_CONF{i}"
                        else:
                            # Single structure, no changes needed
                            new_title = title
                    
                    new_title_map[id(st)] = new_title
                    print(f"    📝 '{st.title}' -> '{new_title}'")
        
        return new_title_map
    
    def update_titles(self, structures, new_title_map):
        """Update structure titles with new tagged names"""
        updated_count = 0
        
        for st in structures:
            new_title = new_title_map.get(id(st))
            if new_title and new_title != st.title:
                old_title = st.title
                st.title = new_title
                
                # Also set as a custom property for tracking
                st.property['s_user_original_title'] = old_title
                st.property['s_user_tagged_title'] = new_title
                
                updated_count += 1
        
        print(f"✏️ Updated {updated_count} structure titles")
    
    def export_to_mae(self, structures):
        """Export structures to MAE file"""
        mae_path = f"{self.output_basename}.mae"
        
        try:
            with StructureWriter(mae_path) as writer:
                for st in structures:
                    writer.append(st)
            
            print(f"✅ MAE file exported to: {mae_path}")
            return mae_path
            
        except Exception as e:
            print(f"❌ Error writing MAE file: {str(e)}")
            return None
    
    def export_to_sdf(self, structures):
        """Export structures to SDF file"""
        sdf_path = f"{self.output_basename}.sdf"
        
        try:
            with StructureWriter(sdf_path) as writer:
                for st in structures:
                    writer.append(st)
            
            print(f"✅ SDF file exported to: {sdf_path}")
            return sdf_path
            
        except Exception as e:
            print(f"❌ Error writing SDF file: {str(e)}")
            return None
    
    def export_to_csv(self, structures):
        """Export structure data to CSV (requires pandas)"""
        if not PANDAS_AVAILABLE:
            print("⚠️ CSV export skipped (pandas not available)")
            return None
        
        rows = []
        
        for st in structures:
            # Start with title
            row = {'Title': st.title}
            
            # Add all structure properties
            for prop_name, prop_value in st.property.items():
                if prop_name not in row:  # Avoid overwriting Title
                    row[prop_name] = prop_value
            
            rows.append(row)
        
        df = pd.DataFrame(rows)
        
        # Ensure Title is first column
        other_cols = [c for c in df.columns if c != 'Title']
        df = df[['Title'] + other_cols]
        
        csv_path = f"{self.output_basename}.csv"
        df.to_csv(csv_path, index=False)
        print(f"✅ CSV exported to: {csv_path}")
        return csv_path
    
    def generate_summary_report(self, structures, grouped_structures):
        """Generate a summary report of the processing"""
        report_lines = []
        report_lines.append("=" * 60)
        report_lines.append("STEREOCHEMISTRY TAGGING SUMMARY REPORT")
        report_lines.append("=" * 60)
        report_lines.append(f"Input file: {self.input_file}")
        report_lines.append(f"Total structures processed: {len(structures)}")
        report_lines.append(f"Unique compound groups: {len(grouped_structures)}")
        report_lines.append("")
        
        # Group analysis
        for title, group in grouped_structures.items():
            if len(group) > 1:
                report_lines.append(f"Group '{title}': {len(group)} structures")
                
                # Show stereochemistry diversity
                stereo_signatures = set()
                for st in group:
                    centers = self.compute_stereochemistry_schrodinger(st)
                    if centers:
                        sig = "_".join([f"{idx}{conf}" for idx, conf in centers])
                        stereo_signatures.add(sig)
                
                if stereo_signatures:
                    report_lines.append(f"  Stereoisomers found: {', '.join(sorted(stereo_signatures))}")
                else:
                    report_lines.append(f"  No stereochemistry detected")
                report_lines.append("")
        
        # Write report
        report_path = f"{self.output_basename}_report.txt"
        with open(report_path, 'w') as f:
            f.write('\n'.join(report_lines))
        
        print(f"📊 Summary report written to: {report_path}")
        return report_path
    
    def run(self, output_format='mae', include_csv=False, include_report=False):
        """Run the complete preprocessing workflow"""
        print(f"🚀 Starting MAE/SDF preprocessing workflow...")
        print(f"📁 Input: {self.input_file}")
        print(f"📝 Output basename: {self.output_basename}")
        print(f"🎯 Output format: {output_format}")
        
        # Step 1: Read structures
        structures = self.read_structures()
        if not structures:
            print("❌ No structures found in input file")
            return None
        
        # Step 2: Normalize titles
        self.normalize_titles(structures)
        
        # Step 3: Group by title
        grouped = self.group_by_title(structures)
        print(f"🔍 Found {len(grouped)} unique compound groups")
        
        # Step 4: Generate new titles with stereochemistry
        new_title_map = self.generate_new_titles(grouped)
        
        # Step 5: Update structure titles
        self.update_titles(structures, new_title_map)
        
        # Step 6: Export results
        results = {}
        
        if output_format.lower() == 'mae':
            results['mae'] = self.export_to_mae(structures)
        elif output_format.lower() == 'sdf':
            results['sdf'] = self.export_to_sdf(structures)
        else:
            # Export both
            results['mae'] = self.export_to_mae(structures)
            results['sdf'] = self.export_to_sdf(structures)
        
        if include_csv:
            results['csv'] = self.export_to_csv(structures)
        
        if include_report:
            results['report'] = self.generate_summary_report(structures, grouped)
        
        return results


def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(
        description="MAE/SDF Preprocessor with Stereochemistry & Conformer Tagging using Schrödinger API",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process MAE file and output MAE
  $SCHRODINGER/run python mae_preprocessor.py input.mae -o tagged_output

  # Process SDF file and output MAE  
  $SCHRODINGER/run python mae_preprocessor.py input.sdf -o output --mae-output
  
  # Include CSV and summary report
  $SCHRODINGER/run python mae_preprocessor.py input.mae -o results --csv --report
        """
    )
    
    parser.add_argument("input_file", help="Input file (.mae, .sdf, .sdf.gz, or .sdfgz)")
    parser.add_argument("-o", "--output", required=True, 
                       help="Output basename (without extension)")
    parser.add_argument("--mae-output", action="store_true",
                       help="Output MAE format (default)")
    parser.add_argument("--sdf-output", action="store_true", 
                       help="Output SDF format")
    parser.add_argument("--csv", action="store_true",
                       help="Also export CSV file with metadata")
    parser.add_argument("--report", action="store_true",
                       help="Generate summary report")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"❌ Error: Input file not found: {args.input_file}")
        sys.exit(1)
    
    # Determine output format
    if args.sdf_output and not args.mae_output:
        output_format = 'sdf'
    elif args.mae_output and not args.sdf_output:
        output_format = 'mae'
    else:
        # Default to MAE
        output_format = 'mae'
    
    # Create preprocessor
    preprocessor = MAEPreprocessor(
        input_file=args.input_file,
        output_basename=args.output
    )
    
    try:
        # Run workflow
        results = preprocessor.run(
            output_format=output_format,
            include_csv=args.csv,
            include_report=args.report
        )
        
        if results:
            print(f"🎉 Workflow completed successfully!")
            for file_type, path in results.items():
                if path:
                    print(f"📁 {file_type.upper()}: {path}")
        else:
            print(f"❌ Workflow failed")
            sys.exit(1)
            
    except Exception as e:
        print(f"❌ Workflow failed: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

import streamlit as st
# Set page configuration
st.set_page_config(
    page_title="GenoScope",
    layout="wide",
    initial_sidebar_state="expanded",
    page_icon="🧬"
)

import os
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from io import BytesIO, StringIO
import subprocess
import tempfile
import base64
import gzip
import json
from collections import defaultdict
import streamlit.components.v1 as components
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# Constants
TEMP_DIR = "temp_results"
os.makedirs(TEMP_DIR, exist_ok=True)

# Custom CSS for better styling
st.markdown("""
<style>
    .main {padding: 2rem;}
    .stButton>button {width: 100%; background-color: #4CAF50; color: white;}
    .stDownloadButton>button {width: 100%; background-color: #2196F3; color: white;}
    .report-view-container {margin-top: -5em;}
    .stProgress > div > div > div > div {background-color: #1f77b4;}
    .st-bb {background-color: transparent;}
    .st-at {background-color: #0d1117;}
    footer {visibility: hidden;}
    .stAlert {padding: 0.5rem;}
    .tool-card {
        border-radius: 10px; 
        padding: 15px; 
        margin-bottom: 15px; 
        box-shadow: 0 4px 8px 0 rgba(0,0,0,0.1);
        background-color: #f8f9fa;
    }
    .tool-card h3 {color: #2c3e50; margin-top: 0;}
    .sidebar .sidebar-content {
        background-color: #f8f9fa;
    }
    h1 {color: #2c3e50;}
    h2 {color: #3498db;}
    .st-bq {border-color: #3498db;}
    .st-cb {background-color: #3498db;}
</style>
""", unsafe_allow_html=True)

# Tool configurations
TOOLS = {
    'Trimming': {
        'Trimmomatic': {
            'command': 'trimmomatic PE -phred33 {input1} {input2} {output1} {output1_unpaired} {output2} {output2_unpaired}',
            'params': {
                'ILLUMINACLIP': 'TruSeq3-PE.fa:2:30:10',
                'LEADING': 3,
                'TRAILING': 3,
                'SLIDINGWINDOW': '4:15',
                'MINLEN': 36
            }
        },
        'Cutadapt': {
            'command': 'cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o {output1} -p {output2} {input1} {input2}',
            'params': {
                'quality-cutoff': 20,
                'minimum-length': 20,
                'overlap': 3
            }
        },
        'fastp': {
            'command': 'fastp -i {input1} -I {input2} -o {output1} -O {output2}',
            'params': {
                'qualified_quality_phred': 15,
                'length_required': 30,
                'cut_front': True,
                'cut_tail': True
            }
        }
    },
    'Alignment': {
        'BWA-MEM': {
            'command': 'bwa mem -t 4 -R "@RG\\tID:{sample}\\tSM:{sample}" {reference} {input1} {input2} > {output}',
            'params': {
                'mark_short_splits': True,
                'min_seed_length': 19,
                'band_width': 100
            }
        },
        'Bowtie2': {
            'command': 'bowtie2 -x {index} -1 {input1} -2 {input2} -S {output}',
            'params': {
                'very-sensitive': True,
                'no-mixed': True,
                'no-discordant': True
            }
        },
        'Minimap2': {
            'command': 'minimap2 -ax sr {reference} {input1} {input2} > {output}',
            'params': {
                'secondary': 'no',
                'min_occurence': 5
            }
        }
    },
    'Variant Calling': {
        'GATK': {
            'command': 'gatk HaplotypeCaller -R {reference} -I {input} -O {output}',
            'params': {
                'stand-call-conf': 20,
                'emit-ref-confidence': 'GVCF'
            }
        },
        'BCFtools': {
            'command': 'bcftools mpileup -Ou -f {reference} {input} | bcftools call -mv -Ob -o {output}',
            'params': {
                'min-MQ': 20,
                'min-BQ': 20
            }
        },
        'FreeBayes': {
            'command': 'freebayes -f {reference} {input} > {output}',
            'params': {
                'min-base-quality': 20,
                'min-alternate-fraction': 0.2
            }
        }
    }
}

# Utility functions
def save_uploaded_file(uploaded_file, directory=TEMP_DIR):
    """Save uploaded file to temporary directory"""
    file_path = os.path.join(directory, uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return file_path

def generate_dummy_fastqc_report(sample_name, read_type):
    """Generate a dummy FastQC report for demonstration"""
    report_path = os.path.join(TEMP_DIR, f"{sample_name}_{read_type}_fastqc.html")
    with open(report_path, "w") as f:
        f.write(f"""
        <html><body>
            <h1>FastQC Report for {sample_name} {read_type}</h1>
            <h2>Basic Statistics</h2>
            <table>
                <tr><th>Measure</th><th>Value</th></tr>
                <tr><td>Total Sequences</td><td>1,500,000</td></tr>
                <tr><td>Sequence Length</td><td>150</td></tr>
                <tr><td>%GC</td><td>52</td></tr>
            </table>
            <h2>Per Base Sequence Quality</h2>
            <img src="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc.png" width="500">
        </body></html>
        """)
    return report_path

def generate_dummy_fastqc_report_2(sample_name, read_type):
    """Generate a dummy FastQC report for demonstration"""
    report_path = os.path.join(TEMP_DIR, f"{sample_name}_{read_type}_fastqc.html")
    with open(report_path, "w") as f:
        f.write(f"""
        <html><body>
            <h1>FastQC Report for {sample_name} {read_type}</h1>
            <h2>Basic Statistics</h2>
            <table>
                <tr><th>Measure</th><th>Value</th></tr>
                <tr><td>Total Sequences</td><td>1,500,000</td></tr>
                <tr><td>Sequence Length</td><td>150</td></tr>
                <tr><td>%GC</td><td>52</td></tr>
            </table>
            <h2>Per Base Sequence Quality</h2>
            <img src="https://i.ytimg.com/vi/YcijCIWuT80/maxresdefault.jpg" width="500">
        </body></html>
        """)
    return report_path

def show_igv(bam_path=None, ref_path=None, vcf_path=None):
    """Display IGV.js genome browser with proper configuration"""
    # Default tracks if no files provided
    if not bam_path:
        bam_path = "https://1000genomes.s3.amazonaws.com/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
        bai_path = "https://1000genomes.s3.amazonaws.com/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai"
    else:
        bai_path = f"{bam_path}.bai"
    
    if not ref_path:
        ref_path = "hg38"
    
    tracks = [
        {
            "name": "Alignments",
            "url": bam_path,
            "indexURL": bai_path,
            "type": "alignment",
            "displayMode": "EXPANDED"
        }
    ]
    
    if vcf_path:
        tracks.append({
            "name": "Variants",
            "url": vcf_path,
            "type": "variant",
            "displayMode": "EXPANDED"
        })
    
    igv_options = {
        "genome": ref_path,
        "tracks": tracks
    }
    
    # IGV HTML component
    igv_html = f"""
    <div id="igv-div" style="padding-top: 20px; padding-bottom: 20px; height: 600px; width: 100%;"></div>
    <script src="https://cdn.jsdelivr.net/npm/igv@2.12.3/dist/igv.min.js"></script>
    <script>
        var options = {json.dumps(igv_options)};
        document.addEventListener("DOMContentLoaded", function() {{
            var igvDiv = document.getElementById("igv-div");
            igv.createBrowser(igvDiv, options)
                .then(function(browser) {{
                    console.log("IGV browser created");
                }});
        }});
    </script>
    """
    
    components.html(igv_html, height=650)

# Pages
def about_page():
    st.title("🧬 Welcome to GenoScope")
    st.image("https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcTTiY1ESXJEr97ZS1aCkB1AkeMhp2mq3IUFdzG29zHIycKHU4bgmCX4LMpsOL8A-gj8WyQ&usqp=CAU", 
             width=800, caption="Next-Generation Sequencing Analysis Platform")
    
    st.markdown("""
    ## Comprehensive Genomic Analysis Platform
    
    GenoScope provides an integrated environment for processing and analyzing next-generation sequencing data.
    
    ### *Key Features:*
    - **End-to-End Pipeline**: FASTQ to VCF with customizable steps
    - **Interactive Visualization**: Explore alignments and variants
    - **Comprehensive Reports**: Quality metrics, variant statistics, and annotations
    - **User-Friendly Interface**: No command-line expertise required

    ### *How to use?*
    - Go to Genomic Pipeline page. Upload the Fastq and reference fasta files. Please make sure the file sizes are less than 200mb.
    - See the fastqc report generated after the files are uploaded.
    - According to the quality and fastqc report, go ahead with trimming, indexing, bam processing, sam processing and variant calling.
    - There are options for tools to be selected and parameters provided which can be customised as per the user's requirement to run the pipeline.
    - This will generate a variant analysis report easy for user to interpret and process.
    - Additionally, a feature to visualise bam file generated against the reference genome using igv has been added which can be accessed through navigation window.
    
    ### *Supported Analyses:*
    - Whole Genome Sequencing
    - Exome Sequencing
    - Targeted Panels
    - RNA-Seq (coming soon)
    
    ### *Contact:*
    - *Name of developer*: Shreya Sancheti
    - *Email*: shreyasancheti190503@gmail.com
    
    ### Acknowledgement:
    - Dr.Kushgra Kashyap, Professor, DESPU
    - LinkedIN ID: https://www.linkedin.com/in/dr-kushagra-kashyap
    """)

def genomic_pipeline_page():
    st.title("🧬 Genomic Analysis Pipeline")
    st.markdown("Configure and execute your genomic analysis workflow from FASTQ to variant calling")
    
    # Initialize session state for pipeline
    if 'pipeline_results' not in st.session_state:
        st.session_state.pipeline_results = None
    if 'current_step' not in st.session_state:
        st.session_state.current_step = "Upload"
    if 'file_paths' not in st.session_state:
        st.session_state.file_paths = {}
    if 'tool_params' not in st.session_state:
        st.session_state.tool_params = {}
    if 'variant_filters' not in st.session_state:
        st.session_state.variant_filters = {
            'min_qual': 30,
            'min_depth': 10,
            'max_af': 0.05,
            'impact': ['HIGH', 'MODERATE'],
            'clin_sig': ['Pathogenic', 'Likely pathogenic']
        }
    
    # Pipeline steps
    steps = ["Upload", "QC & Trimming", "Alignment", "Variant Calling", "Results"]
    current_step = st.session_state.current_step
    
    # Show pipeline progress
    st.subheader("Pipeline Progress")
    cols = st.columns(len(steps))
    for i, step in enumerate(steps):
        with cols[i]:
            if step == current_step:
                st.markdown(f"<div style='background-color: #3498db; color: white; padding: 10px; border-radius: 5px; text-align: center;'><b>{step}</b></div>", unsafe_allow_html=True)
            elif steps.index(current_step) > i:
                st.markdown(f"<div style='background-color: #2ecc71; color: white; padding: 10px; border-radius: 5px; text-align: center;'><b>✓ {step}</b></div>", unsafe_allow_html=True)
            else:
                st.markdown(f"<div style='background-color: #ecf0f1; padding: 10px; border-radius: 5px; text-align: center;'><b>{step}</b></div>", unsafe_allow_html=True)
    
    # Step 1: File Upload
    if current_step == "Upload":
        st.subheader("1. Upload Input Files")
        
        with st.expander("File Upload", expanded=True):
            col1, col2, col3 = st.columns(3)
            with col1:
                r1 = st.file_uploader("Read 1 (FASTQ)", type=['fastq', 'fq', 'fastq.gz', 'fq.gz'])
            with col2:
                r2 = st.file_uploader("Read 2 (FASTQ)", type=['fastq', 'fq', 'fastq.gz', 'fq.gz'])
            with col3:
                ref = st.file_uploader("Reference Genome (FASTA)", type=['fasta', 'fa', 'fasta.gz'])
            
            st.subheader("Sample Information")
            sample_name = st.text_input("Sample Name", "Sample_01")
            st.session_state.sample_name = sample_name
            
            if r1 and r2 and ref and sample_name:
                if st.button("Process Uploads", type="primary"):
                    # Save files
                    st.session_state.file_paths['r1'] = save_uploaded_file(r1)
                    st.session_state.file_paths['r2'] = save_uploaded_file(r2)
                    st.session_state.file_paths['ref'] = save_uploaded_file(ref)
                    
                    # Generate dummy FastQC reports
                    st.session_state.file_paths['fastqc_r1'] = generate_dummy_fastqc_report(sample_name, "R1")
                    st.session_state.file_paths['fastqc_r2'] = generate_dummy_fastqc_report_2(sample_name, "R2")
                    
                    st.session_state.current_step = "QC & Trimming"
                    st.rerun()
    
    # Step 2: QC and Trimming
    elif current_step == "QC & Trimming":
        st.subheader("2. Quality Control & Trimming")
        
        # Show FastQC reports
        with st.expander("Quality Control Reports", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Read 1 QC Report**")
                if os.path.exists(st.session_state.file_paths['fastqc_r1']):
                    with open(st.session_state.file_paths['fastqc_r1'], 'r') as f:
                        components.html(f.read(), height=600, scrolling=True)
            
            with col2:
                st.markdown("**Read 2 QC Report**")
                if os.path.exists(st.session_state.file_paths['fastqc_r2']):
                    with open(st.session_state.file_paths['fastqc_r2'], 'r') as f:
                        components.html(f.read(), height=600, scrolling=True)
        
        # Trimming options
        with st.expander("Trimming Configuration", expanded=True):
            trim_tool = st.selectbox("Select Trimming Tool", list(TOOLS['Trimming'].keys()))
            
            # Show tool-specific parameters
            st.markdown("**Tool Parameters**")
            for param, default in TOOLS['Trimming'][trim_tool]['params'].items():
                if isinstance(default, bool):
                    st.session_state.tool_params[f"trim_{param}"] = st.checkbox(param, default)
                elif isinstance(default, int):
                    st.session_state.tool_params[f"trim_{param}"] = st.number_input(param, value=default)
                else:
                    st.session_state.tool_params[f"trim_{param}"] = st.text_input(param, value=str(default))
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Run Trimming", type="primary"):
                    with st.spinner("Running trimming..."):
                        # In a real app, you would call the actual tool here
                        st.session_state.file_paths['trimmed_r1'] = os.path.join(TEMP_DIR, f"{st.session_state.sample_name}_trimmed_R1.fastq")
                        st.session_state.file_paths['trimmed_r2'] = os.path.join(TEMP_DIR, f"{st.session_state.sample_name}_trimmed_R2.fastq")
                        
                        # Create dummy trimmed files
                        with open(st.session_state.file_paths['trimmed_r1'], 'w') as f:
                            f.write("Dummy trimmed R1 file")
                        with open(st.session_state.file_paths['trimmed_r2'], 'w') as f:
                            f.write("Dummy trimmed R2 file")
                        
                        st.success("Trimming completed!")
                        st.session_state.current_step = "Alignment"
                        st.rerun()
            with col2:
                if st.button("Back to Upload"):
                    st.session_state.current_step = "Upload"
                    st.rerun()
    
    # Step 3: Alignment
    elif current_step == "Alignment":
        st.subheader("3. Alignment")
        
        with st.expander("Alignment Configuration", expanded=True):
            align_tool = st.selectbox("Select Alignment Tool", list(TOOLS['Alignment'].keys()))
            
            # Show tool-specific parameters
            st.markdown("**Alignment Parameters**")
            for param, default in TOOLS['Alignment'][align_tool]['params'].items():
                if isinstance(default, bool):
                    st.session_state.tool_params[f"align_{param}"] = st.checkbox(param, default)
                elif isinstance(default, int):
                    st.session_state.tool_params[f"align_{param}"] = st.number_input(param, value=default)
                else:
                    st.session_state.tool_params[f"align_{param}"] = st.text_input(param, value=str(default))
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Run Alignment", type="primary"):
                    with st.spinner("Running alignment..."):
                        # In a real app, you would call the actual tool here
                        st.session_state.file_paths['aligned_bam'] = os.path.join(TEMP_DIR, f"{st.session_state.sample_name}_aligned.bam")
                        
                        # Create dummy alignment stats
                        st.session_state.pipeline_results = {
                            'alignment_stats': {
                                'total_reads': 15000000,
                                'mapped_reads': 14250000,
                                'mapping_rate': 95.0,
                                'mean_coverage': 45.2,
                                'coverage_uniformity': 0.88
                            }
                        }
                        
                        st.success("Alignment completed!")
                        st.session_state.current_step = "Variant Calling"
                        st.rerun()
            with col2:
                if st.button("Back to Trimming"):
                    st.session_state.current_step = "QC & Trimming"
                    st.rerun()
    
    # Step 4: Variant Calling
    elif current_step == "Variant Calling":
        st.subheader("4. Variant Calling")
        
        with st.expander("Alignment Statistics", expanded=True):
            # Show alignment stats if available
            if st.session_state.pipeline_results and 'alignment_stats' in st.session_state.pipeline_results:
                stats = st.session_state.pipeline_results['alignment_stats']
                col1, col2, col3 = st.columns(3)
                col1.metric("Total Reads", f"{stats['total_reads']:,}")
                col2.metric("Mapped Reads", f"{stats['mapped_reads']:,}", f"{stats['mapping_rate']:.1f}%")
                col3.metric("Mean Coverage", f"{stats['mean_coverage']:.1f}X")
        
        with st.expander("Variant Calling Configuration", expanded=True):
            vc_tool = st.selectbox("Select Variant Caller", list(TOOLS['Variant Calling'].keys()))
            
            # Show tool-specific parameters
            st.markdown("**Variant Calling Parameters**")
            for param, default in TOOLS['Variant Calling'][vc_tool]['params'].items():
                if isinstance(default, bool):
                    st.session_state.tool_params[f"vc_{param}"] = st.checkbox(param, default)
                elif isinstance(default, int):
                    st.session_state.tool_params[f"vc_{param}"] = st.number_input(param, value=default)
                else:
                    st.session_state.tool_params[f"vc_{param}"] = st.text_input(param, value=str(default))
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Run Variant Calling", type="primary"):
                    with st.spinner("Running variant calling..."):
                        # In a real app, you would call the actual tool here
                        st.session_state.file_paths['variants'] = os.path.join(TEMP_DIR, f"{st.session_state.sample_name}_variants.vcf")
                        
                        # Create dummy variant results
                        variants = {
                            'SNP': 1254,
                            'INDEL': 312,
                            'INSERTION': 78,
                            'DELETION': 234,
                            'COMPLEX': 45
                        }
                        
                        annotated_variants = [
                            {'CHROM': '1', 'POS': 12345, 'REF': 'A', 'ALT': 'G', 'TYPE': 'SNP', 
                             'GENE': 'BRCA1', 'EFFECT': 'missense', 'IMPACT': 'HIGH',
                             'CLIN_SIG': 'Pathogenic', 'AF': 0.0012, 'DP': 45, 'QUAL': 120},
                            {'CHROM': '7', 'POS': 552490, 'REF': 'C', 'ALT': 'T', 'TYPE': 'SNP',
                             'GENE': 'CFTR', 'EFFECT': 'synonymous', 'IMPACT': 'LOW',
                             'CLIN_SIG': 'Benign', 'AF': 0.42, 'DP': 38, 'QUAL': 110},
                            {'CHROM': '17', 'POS': 412760, 'REF': 'G', 'ALT': 'GA', 'TYPE': 'INDEL',
                             'GENE': 'TP53', 'EFFECT': 'frameshift', 'IMPACT': 'HIGH',
                             'CLIN_SIG': 'Pathogenic', 'AF': 0.0001, 'DP': 52, 'QUAL': 150},
                            {'CHROM': '13', 'POS': 329140, 'REF': 'C', 'ALT': 'T', 'TYPE': 'SNP',
                             'GENE': 'BRCA2', 'EFFECT': 'missense', 'IMPACT': 'MODERATE',
                             'CLIN_SIG': 'Likely pathogenic', 'AF': 0.0008, 'DP': 48, 'QUAL': 135},
                            {'CHROM': '9', 'POS': 107620, 'REF': 'G', 'ALT': 'A', 'TYPE': 'SNP',
                             'GENE': 'FANCC', 'EFFECT': 'synonymous', 'IMPACT': 'LOW',
                             'CLIN_SIG': 'Benign', 'AF': 0.35, 'DP': 42, 'QUAL': 105}
                        ]
                        
                        st.session_state.pipeline_results.update({
                            'variant_stats': variants,
                            'annotated_variants': annotated_variants,
                            'all_variants': pd.DataFrame(annotated_variants)
                        })
                        
                        st.success("Variant calling completed!")
                        st.session_state.current_step = "Results"
                        st.rerun()
            with col2:
                if st.button("Back to Alignment"):
                    st.session_state.current_step = "Alignment"
                    st.rerun()
    
    # Step 5: Results
    elif current_step == "Results":
        st.subheader("5. Analysis Results")
        
        if st.session_state.pipeline_results:
            # Variant statistics
            st.markdown("### Variant Statistics")
            variants = st.session_state.pipeline_results.get('variant_stats', {})
            total_variants = sum(variants.values())
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Variants", f"{total_variants:,}", delta_color="off")
            col2.metric("SNPs", f"{variants.get('SNP', 0):,}", 
                       f"{variants.get('SNP', 0)/total_variants*100:.1f}%" if total_variants else "0%")
            col3.metric("Indels", f"{variants.get('INDEL', 0) + variants.get('INSERTION', 0) + variants.get('DELETION', 0):,}",
                       f"{(variants.get('INDEL', 0) + variants.get('INSERTION', 0) + variants.get('DELETION', 0))/total_variants*100:.1f}%" if total_variants else "0%")
            
            # Variant type pie chart
            variant_types = ['SNP', 'INDEL', 'INSERTION', 'DELETION', 'COMPLEX']
            counts = [variants.get(t, 0) for t in variant_types]
            
            fig = px.pie(names=variant_types, values=counts,
                        title="Variant Types", hole=0.4,
                        color_discrete_sequence=px.colors.qualitative.Pastel)
            st.plotly_chart(fig, use_container_width=True)
            
            # Variant filtering
            with st.expander("Variant Filters", expanded=True):
                st.session_state.variant_filters['min_qual'] = st.slider("Minimum Quality", 0, 200, st.session_state.variant_filters['min_qual'])
                st.session_state.variant_filters['min_depth'] = st.slider("Minimum Depth", 0, 100, st.session_state.variant_filters['min_depth'])
                st.session_state.variant_filters['max_af'] = st.slider("Maximum Allele Frequency", 0.0, 1.0, st.session_state.variant_filters['max_af'])
                st.session_state.variant_filters['impact'] = st.multiselect("Impact", ['HIGH', 'MODERATE', 'LOW', 'MODIFIER'], st.session_state.variant_filters['impact'])
                st.session_state.variant_filters['clin_sig'] = st.multiselect("Clinical Significance", 
                                             ['Pathogenic', 'Likely pathogenic', 'Uncertain significance', 
                                              'Likely benign', 'Benign', 'Conflicting'],
                                             st.session_state.variant_filters['clin_sig'])
            
            # Apply filters
            filtered_variants = st.session_state.pipeline_results['all_variants'][
                (st.session_state.pipeline_results['all_variants']['QUAL'] >= st.session_state.variant_filters['min_qual']) &
                (st.session_state.pipeline_results['all_variants']['DP'] >= st.session_state.variant_filters['min_depth']) &
                (st.session_state.pipeline_results['all_variants']['AF'] <= st.session_state.variant_filters['max_af']) &
                (st.session_state.pipeline_results['all_variants']['IMPACT'].isin(st.session_state.variant_filters['impact'])) &
                (st.session_state.pipeline_results['all_variants']['CLIN_SIG'].isin(st.session_state.variant_filters['clin_sig']))
            ]
            
            # Display filtered variants
            st.markdown(f"### Filtered Variants ({len(filtered_variants)} found)")
            st.dataframe(filtered_variants)
            
            # Variant visualizations
            st.markdown("### Variant Visualizations")
            
            tab1, tab2, tab3 = st.tabs(["By Chromosome", "By Impact", "Clinical Significance"])
            
            with tab1:
                chrom_counts = filtered_variants['CHROM'].value_counts().reset_index()
                chrom_counts.columns = ['Chromosome', 'Count']
                fig = px.bar(chrom_counts, x='Chromosome', y='Count', 
                            title="Variants by Chromosome",
                            color='Chromosome',
                            color_discrete_sequence=px.colors.qualitative.Pastel)
                st.plotly_chart(fig, use_container_width=True)
            
            with tab2:
                impact_counts = filtered_variants['IMPACT'].value_counts().reset_index()
                impact_counts.columns = ['Impact', 'Count']
                fig = px.pie(impact_counts, names='Impact', values='Count',
                            title="Variants by Impact",
                            color='Impact',
                            color_discrete_map={
                                'HIGH': '#FF7F0E',
                                'MODERATE': '#1F77B4',
                                'LOW': '#2CA02C',
                                'MODIFIER': '#D62728'
                            })
                st.plotly_chart(fig, use_container_width=True)
            
            with tab3:
                clin_counts = filtered_variants['CLIN_SIG'].value_counts().reset_index()
                clin_counts.columns = ['Clinical Significance', 'Count']
                fig = px.bar(clin_counts, x='Clinical Significance', y='Count',
                            color='Clinical Significance',
                            title="Clinical Significance of Variants",
                            color_discrete_map={
                                'Pathogenic': '#D62728',
                                'Likely pathogenic': '#FF7F0E',
                                'Uncertain significance': '#8C564B',
                                'Likely benign': '#1F77B4',
                                'Benign': '#2CA02C'
                            })
                st.plotly_chart(fig, use_container_width=True)
            
            # Visualization
            st.markdown("### Genome Browser")
            if st.button("View in IGV", type="primary"):
                show_igv(
                    bam_path=st.session_state.file_paths.get('aligned_bam', ""),
                    ref_path=st.session_state.file_paths.get('ref', ""),
                    vcf_path=st.session_state.file_paths.get('variants', "")
                )
            
            # Download results
            st.markdown("### Download Results")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.download_button(
                    label="Download Variants (VCF)",
                    data="Dummy VCF data would be here in a real implementation",
                    file_name=f"{st.session_state.sample_name}_variants.vcf",
                    mime="text/vcf"
                )
            with col2:
                st.download_button(
                    label="Download Alignment Stats",
                    data=json.dumps(st.session_state.pipeline_results.get('alignment_stats', {})),
                    file_name=f"{st.session_state.sample_name}_alignment_stats.json",
                    mime="application/json"
                )
            with col3:
                st.download_button(
                    label="Download Full Report",
                    data="Dummy full report data would be here",
                    file_name=f"{st.session_state.sample_name}_report.pdf",
                    mime="application/pdf"
                )
        
        if st.button("Restart Pipeline", type="primary"):
            st.session_state.current_step = "Upload"
            st.rerun()

def visualization_page():
    st.title("🔬 Genomic Visualization")
    st.markdown("Interactive visualization of genomic data")
    
    # Visualization options
    viz_option = st.selectbox("Select Visualization Type", 
                             ["IGV Browser", "Coverage Plot", "Variant Distribution"])
    
    if viz_option == "IGV Browser":
        st.markdown("### Interactive Genome Browser")
        
        # For a real implementation, you would use actual BAM/VCF files
        st.info("In a real implementation, this would show actual BAM/VCF data")
        
        if st.button("Launch IGV with Demo Data", type="primary"):
            show_igv()
    
    elif viz_option == "Coverage Plot":
        st.markdown("### Coverage Across Genomic Regions")
        
        # Generate dummy coverage data
        positions = list(range(1, 1001))
        coverage = [max(0, int(50 + 30 * np.sin(x/50) + np.random.normal(0, 5))) for x in positions]
        
        fig = px.line(x=positions, y=coverage, 
                     title="Coverage Across Genomic Region",
                     labels={'x': 'Genomic Position', 'y': 'Coverage Depth'},
                     color_discrete_sequence=['#3498db'])
        st.plotly_chart(fig, use_container_width=True)
    
    elif viz_option == "Variant Distribution":
        st.markdown("### Variant Distribution Across Genome")
        
        # Generate dummy variant data
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']
        variant_counts = {chr: np.random.randint(5, 100) for chr in chromosomes}
        variant_counts['chr1'] += 150  # Make chr1 more prominent
        
        fig = px.bar(x=list(variant_counts.keys()), y=list(variant_counts.values()),
                    title="Variants by Chromosome",
                    labels={'x': 'Chromosome', 'y': 'Variant Count'},
                    color=list(variant_counts.keys()),
                    color_discrete_sequence=px.colors.qualitative.Pastel)
        fig.update_layout(xaxis_categoryorder='total descending', showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

# Main app
def main():
    # Sidebar navigation
    st.sidebar.title("GenoScope")
    st.sidebar.image("https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcTTiY1ESXJEr97ZS1aCkB1AkeMhp2mq3IUFdzG29zHIycKHU4bgmCX4LMpsOL8A-gj8WyQ&usqp=CAU", 
                    use_column_width=True)
    
    pages = {
        "About": about_page,
        "Genomic Pipeline": genomic_pipeline_page,
        "Visualization": visualization_page
    }
    
    selection = st.sidebar.radio("Navigation", list(pages.keys()))
    pages[selection]()

if __name__ == "__main__":
    main()

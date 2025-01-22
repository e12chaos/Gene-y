import streamlit as st
import subprocess
import os
import glob
import pandas as pd
from pathlib import Path

# --------------------------- Session State Initialization --------------------------- #
# Initialize session state for pipeline status
if 'pipeline_status' not in st.session_state:
    st.session_state['pipeline_status'] = {
        'fastqc': False,
        'trimming': False,
        'alignment': False,
        'variant_calling': False
    }

# Initialize session state for pipeline settings
if 'settings' not in st.session_state:
    st.session_state['settings'] = {
        'run_qc': True,
        'run_trim': True,
        'run_align': True,
        'run_variant': True,
        'threads': 4
    }

# Initialize session state for file paths
if 'files' not in st.session_state:
    st.session_state['files'] = {
        'reference_genome': "data/reference/hg38.fa",
        'adapter_file': None,
        'fastq_files': []
    }

# --------------------------- Pipeline Runner Class --------------------------- #
class PipelineRunner:
    def __init__(self):
        self.setup_directories()

    def setup_directories(self):
        """Create necessary directories for the pipeline."""
        directories = [
            "data/raw", "data/reference", "data/adapters",
            "output/QC_Reports", "output/trim",
            "output/alignment", "output/variants",
            "output/VC/VC_sorting", "output/VC/VC_dedup_recal",
            "output/VC/VC_dedup_recal/dedup_recal_metrics",
            "output/VC/VC_bqsr_call/vcf",
            "logs"
        ]
        for dir_path in directories:
            os.makedirs(dir_path, exist_ok=True)

    def run_command(self, command, progress_bar=None):
        """
        Execute a shell command and handle its output.

        Args:
            command (str): The command to execute.
            progress_bar (st.progress): Streamlit progress bar.

        Returns:
            tuple: (stdout, stderr)
        """
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            universal_newlines=True
        )

        output = []
        error = []

        while True:
            stdout_line = process.stdout.readline()
            stderr_line = process.stderr.readline()

            if stdout_line == '' and stderr_line == '' and process.poll() is not None:
                break

            if stdout_line:
                output.append(stdout_line)
                if progress_bar:
                    # Update progress bar based on output lines (placeholder logic)
                    progress = min(len(output) / 100, 1.0)
                    progress_bar.progress(progress)

            if stderr_line:
                error.append(stderr_line)
                st.warning(stderr_line.strip())

        return ''.join(output), ''.join(error)

    def save_uploaded_files(self, uploaded_files, directory):
        """
        Save uploaded files to the specified directory.

        Args:
            uploaded_files (UploadedFile or list): Uploaded file(s) from Streamlit.
            directory (str): Target directory to save files.

        Returns:
            str or list: Path(s) to the saved file(s).
        """
        if not isinstance(uploaded_files, list):
            uploaded_files = [uploaded_files]
        saved_paths = []
        for uploaded_file in uploaded_files:
            if uploaded_file is not None:
                save_path = Path(directory) / uploaded_file.name
                with open(save_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                saved_paths.append(str(save_path))
        return saved_paths[0] if len(saved_paths) == 1 else saved_paths

    def run_pipeline_step(self, step, input_files, adapter_file=None, reference_file=None, known_sites_vcf=None):
        """
        Execute a specific pipeline step.

        Args:
            step (str): The pipeline step to execute.
            input_files (list): List of input file paths.
            adapter_file (str, optional): Path to the adapter file.
            reference_file (str, optional): Path to the reference genome.
            known_sites_vcf (str, optional): Path to the known sites VCF file (for variant calling).

        Returns:
            str or None: Output from the command or None if failed.
        """
        script_path = f"src/{step}.sh"

        # For variant_calling, pass sample_id, reference_genome, and known_sites_vcf
        if step == "variant_calling":
            # Assuming input_files are the final BAM files after BQSR
            # Extract sample IDs from filenames
            for bam_path in input_files:
                sample_id = Path(bam_path).stem.replace('_bqsr', '')  # e.g., BM002B_A_S1_L001_bqsr -> BM002B_A_S1_L001
                command = f"bash {script_path} {sample_id} {reference_file} {known_sites_vcf}"
                progress_bar = st.progress(0)
                status_text = st.empty()

                try:
                    status_text.text(f"Running {step.capitalize()} for sample {sample_id}...")
                    output, error = self.run_command(command, progress_bar)

                    if error and "error" in error.lower():
                        raise Exception(error)

                    st.session_state['pipeline_status'][step] = True
                    status_text.text(f"{step.capitalize()} completed successfully for sample {sample_id}!")
                except Exception as e:
                    st.error(f"Error in {step} for sample {sample_id}: {str(e)}")
                    return None
                finally:
                    progress_bar.empty()
        else:
            # For other steps, pass input_files and other necessary arguments
            input_files_str = " ".join(input_files)
            if step == "fastqc":
                command = f"bash {script_path} {input_files_str}"
            elif step == "trimming":
                command = f"bash {script_path} {adapter_file} {input_files_str}"
            elif step == "alignment":
                command = f"bash {script_path} {reference_file} -t {st.session_state['settings']['threads']} {input_files_str}"
            else:
                command = ""

            progress_bar = st.progress(0)
            status_text = st.empty()

            try:
                status_text.text(f"Running {step.capitalize()}...")
                output, error = self.run_command(command, progress_bar)

                if error and "error" in error.lower():
                    raise Exception(error)

                st.session_state['pipeline_status'][step] = True
                status_text.text(f"{step.capitalize()} completed successfully!")
                return output

            except Exception as e:
                st.error(f"Error in {step}: {str(e)}")
                return None
            finally:
                progress_bar.empty()

    def run_full_pipeline(self, fastq_files, reference_file, adapter_file, known_sites_vcf):
        """
        Execute the entire pipeline based on selected settings.

        Args:
            fastq_files (list): List of FASTQ file paths.
            reference_file (str): Path to the reference genome.
            adapter_file (str): Path to the adapter file.
            known_sites_vcf (str): Path to the known sites VCF file.

        Returns:
            bool: True if pipeline completes successfully, else False.
        """
        settings = st.session_state['settings']
        steps = []

        if settings['run_qc']:
            steps.append(('fastqc', fastq_files, None, None, None))
        if settings['run_trim']:
            steps.append(('trimming', fastq_files, adapter_file, None, None))
        if settings['run_align']:
            trimmed_files = glob.glob("output/trim/*_paired.fastq.gz")
            steps.append(('alignment', trimmed_files, None, reference_file, None))
        if settings['run_variant']:
            aligned_bam_files = glob.glob("output/alignment/*_bqsr.bam")  # Assuming alignment step outputs BQSR BAMs
            steps.append(('variant_calling', aligned_bam_files, None, reference_file, known_sites_vcf))

        progress_bar = st.progress(0)
        status_text = st.empty()

        total_steps = len(steps)
        for idx, (step, inputs, adapter, reference, known_vcf) in enumerate(steps):
            status_text.text(f"Running {step.capitalize()} ({idx + 1}/{total_steps})")
            if step == "variant_calling":
                result = self.run_pipeline_step(step, inputs, adapter, reference, known_vcf)
            else:
                result = self.run_pipeline_step(step, inputs, adapter, reference)
            if result is None and step != "variant_calling":
                return False
            progress_bar.progress((idx + 1) / total_steps)

        status_text.text("Pipeline completed successfully!")
        return True

# --------------------------- UI Rendering Functions --------------------------- #
def render_autonomous_mode():
    """Render the Autonomous Pipeline Mode UI."""
    st.header("Autonomous Pipeline Mode")
    st.markdown("### Required Files")

    # Reference Genome Section
    st.markdown("#### Reference Genome")
    
    # Extract the reference genome name without path and extension
    default_ref_path = Path(st.session_state['files']['reference_genome'])
    default_ref_name = default_ref_path.stem  # Extracts 'hg38' from 'hg38.fa'
    
    st.write(f"**Default Reference Genome:** {default_ref_name}")
    
    # Checkbox to use a different reference genome
    use_custom_ref = st.checkbox("Use a different Reference Genome", key="use_custom_ref_autonomous")
    
    if use_custom_ref:
        uploaded_reference = st.file_uploader(
            "Upload a different Reference Genome (optional)", 
            type=["fa", "fasta"],
            help="If you want to use a different reference genome, upload it here. Otherwise, the default (hg38.fa) will be used.",
            key="upload_reference_autonomous"
        )
    else:
        uploaded_reference = None

    # Adapter File Uploader
    adapter = st.file_uploader(
        "Adapter File", 
        type=["fa", "fasta"],
        help="Upload your adapter sequences file.",
        key="adapter_autonomous"
    )

    # FASTQ Files Uploader
    fastq_files = st.file_uploader(
        "FASTQ Files", 
        type=["fastq", "fq", "fastq.gz", "fq.gz"],
        accept_multiple_files=True,
        help="Upload your FASTQ files for analysis.",
        key="fastq_autonomous"
    )
    
    # Known Sites VCF File Uploader
    known_sites_vcf = st.file_uploader(
        "Known Sites VCF File", 
        type=["vcf"],
        help="Upload the known sites VCF file required for variant recalibration.",
        key="known_sites_vcf_autonomous"
    )
    
    # Start Pipeline Button
    if adapter and fastq_files and known_sites_vcf:
        pipeline = PipelineRunner()
        if st.button("Start Pipeline", type="primary", key="start_pipeline_autonomous"):
            # Handle Reference Genome
            if use_custom_ref and uploaded_reference:
                reference_path = pipeline.save_uploaded_files(uploaded_reference, "data/reference")
                st.session_state['files']['reference_genome'] = reference_path
            else:
                reference_path = st.session_state['files']['reference_genome']
            
            # Save Adapter File
            adapter_path = pipeline.save_uploaded_files(adapter, "data/adapters")
            st.session_state['files']['adapter_file'] = adapter_path
            
            # Save FASTQ Files
            fastq_paths = pipeline.save_uploaded_files(fastq_files, "data/raw")
            st.session_state['files']['fastq_files'] = fastq_paths
            
            # Save Known Sites VCF
            known_vcf_path = pipeline.save_uploaded_files(known_sites_vcf, "data/known_sites")
            
            # Run the full pipeline
            pipeline.run_full_pipeline(fastq_paths, reference_path, adapter_path, known_vcf_path)

def render_interactive_mode():
    """Render the Interactive Pipeline Mode UI."""
    st.header("Interactive Pipeline Mode")
    
    selected_step = st.radio(
        "Select Pipeline Step", 
        ["fastqc", "trimming", "alignment", "variant_calling"],
        horizontal=True,
        key="selected_step_interactive"
    )
    
    if selected_step == "fastqc":
        st.markdown("### Upload FASTQ Files for FastQC")
        fastq_files = st.file_uploader(
            "Upload FASTQ files", 
            type=["fastq", "fq", "fastq.gz", "fq.gz"], 
            accept_multiple_files=True,
            help="Upload your FASTQ files for Quality Control.",
            key="fastqc_fastq_interactive"
        )
        if fastq_files and st.button("Run FastQC", key="run_fastqc_interactive"):
            pipeline = PipelineRunner()
            saved_files = pipeline.save_uploaded_files(fastq_files, "data/raw")
            pipeline.run_pipeline_step('fastqc', saved_files)
            
    elif selected_step == "trimming":
        st.markdown("### Upload Adapter and FASTQ Files for Trimming")
        adapter = st.file_uploader(
            "Upload Adapter File", 
            type=["fa", "fasta"], 
            help="Upload your adapter sequences file.",
            key="trimming_adapter_interactive"
        )
        fastq_files = st.file_uploader(
            "Upload FASTQ files", 
            type=["fastq", "fq", "fastq.gz", "fq.gz"], 
            accept_multiple_files=True,
            help="Upload your FASTQ files for trimming.",
            key="trimming_fastq_interactive"
        )
        
        if adapter and fastq_files and st.button("Run Trimming", key="run_trimming_interactive"):
            pipeline = PipelineRunner()
            adapter_path = pipeline.save_uploaded_files(adapter, "data/adapters")
            fastq_paths = pipeline.save_uploaded_files(fastq_files, "data/raw")
            pipeline.run_pipeline_step('trimming', fastq_paths, adapter_file=adapter_path)
            
    elif selected_step == "alignment":
        st.markdown("### Upload Reference Genome and Trimmed FASTQ Files for Alignment")
        st.write(f"**Current Reference Genome:** {Path(st.session_state['files']['reference_genome']).stem}")
        use_custom_ref = st.checkbox("Use a different Reference Genome", key="use_custom_ref_alignment_interactive")
        
        if use_custom_ref:
            uploaded_reference = st.file_uploader(
                "Upload a different Reference Genome (optional)", 
                type=["fa", "fasta"],
                help="If you want to use a different reference genome, upload it here. Otherwise, the current reference will be used.",
                key="upload_reference_alignment_interactive"
            )
        else:
            uploaded_reference = None
        
        fastq_files = st.file_uploader(
            "Upload Trimmed FASTQ files", 
            type=["fastq", "fq", "fastq.gz", "fq.gz"], 
            accept_multiple_files=True,
            help="Upload your trimmed FASTQ files for alignment.",
            key="alignment_fastq_interactive"
        )
        
        if fastq_files and st.button("Run Alignment", key="run_alignment_interactive"):
            pipeline = PipelineRunner()
            # Handle Reference Genome
            if use_custom_ref and uploaded_reference:
                reference_path = pipeline.save_uploaded_files(uploaded_reference, "data/reference")
                st.session_state['files']['reference_genome'] = reference_path
            else:
                reference_path = st.session_state['files']['reference_genome']
            
            # Save Trimmed FASTQ Files
            fastq_paths = pipeline.save_uploaded_files(fastq_files, "data/raw")
            pipeline.run_pipeline_step('alignment', fastq_paths, reference_file=reference_path)
            
    elif selected_step == "variant_calling":
        st.markdown("### Upload BAM/SAM Files and Known Sites VCF for Variant Calling")
        bam_files = st.file_uploader(
            "Upload BAM/SAM files", 
            type=["bam", "sam"], 
            accept_multiple_files=True,
            help="Upload your alignment files for variant calling.",
            key="variant_calling_bam_interactive"
        )
        known_sites_vcf = st.file_uploader(
            "Upload Known Sites VCF File", 
            type=["vcf"],
            help="Upload the known sites VCF file required for variant recalibration.",
            key="known_sites_vcf_interactive"
        )
        
        if bam_files and known_sites_vcf and st.button("Run Variant Calling", key="run_variant_calling_interactive"):
            pipeline = PipelineRunner()
            saved_files = pipeline.save_uploaded_files(bam_files, "data/alignment")
            saved_known_vcf = pipeline.save_uploaded_files(known_sites_vcf, "data/known_sites")
            pipeline.run_pipeline_step('variant_calling', saved_files, reference_file=st.session_state['files']['reference_genome'], known_sites_vcf=saved_known_vcf)

def render_results_section():
    """Render the Results section with tabs for each pipeline step."""
    st.header("Results")
    
    tabs = st.tabs(["🧬 Quality Control", "✂️ Trimming", "🔗 Alignment", "🔬 Variant Calling"])
    
    with tabs[0]:
        st.subheader("Quality Control Results")
        qc_files = glob.glob("output/QC_Reports/*")
        if qc_files:
            for file in qc_files:
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.markdown(f"**{Path(file).name}**")
                with col2:
                    with open(file, "rb") as f:
                        st.download_button(
                            "Download", 
                            f, 
                            file_name=Path(file).name,
                            mime="text/html" if file.endswith('.html') else "text/plain",
                            key=f"download_qc_{Path(file).name}"
                        )
        else:
            st.info("No Quality Control reports available.")
    
    with tabs[1]:
        st.subheader("Trimming Results")
        trim_files = glob.glob("output/trim/*.fastq.gz")
        if trim_files:
            trim_stats = pd.DataFrame({
                'File': [Path(f).name for f in trim_files],
                'Size': [f"{os.path.getsize(f)/1024/1024:.2f} MB" for f in trim_files]
            })
            st.dataframe(trim_stats)
        else:
            st.info("No trimming results available.")
    
    with tabs[2]:
        st.subheader("Alignment Results")
        align_files = glob.glob("output/alignment/*.{sam,bam}")
        if align_files:
            for file in align_files:
                st.markdown(f"**{Path(file).name}**")
                if st.button(f"View Stats - {Path(file).name}", key=f"view_stats_{Path(file).name}"):
                    try:
                        stats = subprocess.check_output(f"samtools flagstat {file}", shell=True).decode()
                        st.code(stats)
                    except subprocess.CalledProcessError as e:
                        st.error(f"Error retrieving stats: {e}")
        else:
            st.info("No alignment results available.")
    
    with tabs[3]:
        st.subheader("Variant Calling Results")
        vcf_files = glob.glob("output/VC/VC_bqsr_call/vcf/*.vcf")
        if vcf_files:
            for file in vcf_files:
                st.markdown(f"**{Path(file).name}**")
                if st.button(f"View Summary - {Path(file).name}", key=f"view_summary_{Path(file).name}"):
                    try:
                        stats = subprocess.check_output(f"bcftools stats {file}", shell=True).decode()
                        st.code(stats)
                    except subprocess.CalledProcessError as e:
                        st.error(f"Error retrieving summary: {e}")
        else:
            st.info("No variant calling results available.")

# --------------------------- Main Function --------------------------- #
def main():
    """Main function to run the Streamlit app."""
    st.set_page_config(
        page_title="Gene-y: Human WGS Pipeline (hg38)",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    st.title("Gene-y: Human WGS Pipeline (hg38)")
    
    # Select Pipeline Mode
    mode = st.radio(
        "Pipeline Mode",
        ["Autonomous", "Interactive"],
        horizontal=True,
        key="pipeline_mode_selection"
    )
    
    # Conditional Sidebar Content
    if mode == "Autonomous":
        with st.sidebar:
            # Centering content in the sidebar using HTML
            st.markdown("<div style='text-align: center;'>", unsafe_allow_html=True)
            st.markdown("### Pipeline Status")
            for step_name, completed in st.session_state['pipeline_status'].items():
                st.markdown(f"{'✅' if completed else '⏳'} {step_name.capitalize()}")
    
            st.markdown("### Pipeline Settings")
            st.session_state['settings']['run_qc'] = st.checkbox("Run FastQC", key="run_qc_autonomous", value=True)
            st.session_state['settings']['run_trim'] = st.checkbox("Run Trimming", key="run_trim_autonomous", value=True)
            st.session_state['settings']['run_align'] = st.checkbox("Run Alignment", key="run_align_autonomous", value=True)
            st.session_state['settings']['run_variant'] = st.checkbox("Run Variant Calling", key="run_variant_autonomous", value=True)
            st.session_state['settings']['threads'] = st.number_input(
                "Threads", 
                min_value=1, 
                max_value=32, 
                value=4,
                step=1,
                key="threads_autonomous"
            )
            st.markdown("</div>", unsafe_allow_html=True)
    else:
        with st.sidebar:
            # Displaying "Interactive Manual Mode" in the sidebar
            st.markdown("<div style='text-align: center;'>", unsafe_allow_html=True)
            st.markdown("### Interactive Manual Mode")
            st.markdown("</div>", unsafe_allow_html=True)

    # Render the selected mode
    if mode == "Autonomous":
        render_autonomous_mode()
    else:
        render_interactive_mode()
        
    # Render the results section
    st.markdown("---")
    render_results_section()

# --------------------------- Entry Point --------------------------- #
if __name__ == "__main__":
    main()

import streamlit as st
import subprocess
import os
import glob

def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    return output, error

def run_pipeline_step(step, input_files, adapter_file=None, reference_file=None):
    script_path = f"src/{step}.sh"
    input_files_str = " ".join(input_files)
    
    if step == "trimming":
        command = f"bash {script_path} {adapter_file} {input_files_str}"
    elif step == "alignment":
        command = f"bash {script_path} {reference_file} {input_files_str}"
    else:
        command = f"bash {script_path} {input_files_str}"
    
    output, error = run_command(command)
    return output, error

st.title("Gene-y: Human WGS Pipeline (hg38)")

# Ensure the necessary directories exist
os.makedirs("data/raw", exist_ok=True)
os.makedirs("output/QC_Reports", exist_ok=True)
os.makedirs("output/trim", exist_ok=True)
os.makedirs("output/alignment", exist_ok=True)

# Sidebar for pipeline steps
st.sidebar.header("Pipeline Steps")
steps = ["fastqc", "trimming", "alignment", "variant_calling"]
selected_step = st.sidebar.radio("Select a Step to Run", steps)

# Main screen content based on selected step
if selected_step == "fastqc":
    st.header("FastQC Step")
    uploaded_files = st.file_uploader("Upload FASTQ files for QC", type="fastq", accept_multiple_files=True)
    
    if uploaded_files:
        saved_files = []
        for uploaded_file in uploaded_files:
            save_path = os.path.join("data/raw", uploaded_file.name)
            with open(save_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            saved_files.append(save_path)
        st.success("Files have been uploaded successfully.")
        
        if st.button("Run FastQC"):
            st.write("Running FastQC step...")
            output, error = run_pipeline_step(selected_step, saved_files)
            if error:
                st.error(f"Error in FastQC step: {error.decode()}")
            else:
                st.success("FastQC step completed successfully.")
                st.text(output.decode())

elif selected_step == "trimming":
    st.header("Trimming Step")
    adapter_file = st.file_uploader("Upload Adapter File", type="fa")
    uploaded_files = st.file_uploader("Upload FASTQ files for Trimming", type="fastq", accept_multiple_files=True)

    if uploaded_files and adapter_file:
        saved_files = []
        for uploaded_file in uploaded_files:
            save_path = os.path.join("data/raw", uploaded_file.name)
            with open(save_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            saved_files.append(save_path)

        adapter_save_path = os.path.join("data", adapter_file.name)
        with open(adapter_save_path, "wb") as f:
            f.write(adapter_file.getbuffer())
        st.success("Adapter file and FASTQ files have been uploaded successfully.")

        if st.button("Run Trimming"):
            st.write("Running trimming step...")
            output, error = run_pipeline_step(selected_step, saved_files, adapter_save_path)
            if error:
                st.error(f"Error in trimming step: {error.decode()}")
            else:
                st.success("Trimming step completed successfully.")
                st.text(output.decode())

elif selected_step == "alignment":
    st.header("Alignment Step")
    reference_file = st.file_uploader("Upload Reference Genome (hg38.fa)", type="fa")
    uploaded_files = st.file_uploader("Upload FASTQ files for Alignment", type="fastq", accept_multiple_files=True)

    if uploaded_files and reference_file:
        saved_files = []
        for uploaded_file in uploaded_files:
            save_path = os.path.join("data/raw", uploaded_file.name)
            with open(save_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            saved_files.append(save_path)

        reference_save_path = os.path.join("data", reference_file.name)
        with open(reference_save_path, "wb") as f:
            f.write(reference_file.getbuffer())
        st.success("Reference genome and FASTQ files have been uploaded successfully.")

        if st.button("Run Alignment"):
            st.write("Running alignment step...")
            output, error = run_pipeline_step(selected_step, saved_files, reference_file=reference_save_path)
            if error:
                st.error(f"Error in alignment step: {error.decode()}")
            else:
                st.success("Alignment step completed successfully.")
                st.text(output.decode())

elif selected_step == "variant_calling":
    st.header("Variant Calling Step")
    uploaded_files = st.file_uploader("Upload SAM/BAM files for Variant Calling", type=["sam", "bam"], accept_multiple_files=True)
    
    if uploaded_files:
        saved_files = []
        for uploaded_file in uploaded_files:
            save_path = os.path.join("data/alignment", uploaded_file.name)
            with open(save_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            saved_files.append(save_path)
        st.success("SAM/BAM files have been uploaded successfully.")

        if st.button("Run Variant Calling"):
            st.write("Running variant calling step...")
            output, error = run_pipeline_step(selected_step, saved_files)
            if error:
                st.error(f"Error in variant calling step: {error.decode()}")
            else:
                st.success("Variant calling step completed successfully.")
                st.text(output.decode())

# Results Tabbed Section
st.header("Results")

# Create tabs for different results categories
tabs = st.tabs(["🧬 Quality Control", "✂️ Trimming", "🔗 Alignment", "🔬 Variant Calling"])

# Tab 1: Quality Control Results
with tabs[0]:
    st.subheader("Quality Control Results")
    qc_files = glob.glob("output/QC_Reports/*.html")
    
    for qc_file in qc_files:
        st.markdown(f"**{os.path.basename(qc_file)}**")
        st.progress(100)  # Example progress bar (replace with actual progress logic if available)
        st.download_button(
            label="Download",
            data=open(qc_file, "rb"),
            file_name=os.path.basename(qc_file),
            mime="text/html"
        )

# Tab 2: Trimming Results
with tabs[1]:
    st.subheader("Trimming Results")
    trimming_files = glob.glob("output/trim/*.fastq.gz")
    
    for trimming_file in trimming_files:
        st.markdown(f"**{os.path.basename(trimming_file)}**")
        st.progress(100)  # Example progress bar (replace with actual progress logic if available)
        st.download_button(
            label="Download",
            data=open(trimming_file, "rb"),
            file_name=os.path.basename(trimming_file),
            mime="application/gzip"
        )

# Tab 3: Alignment Results
with tabs[2]:
    st.subheader("Alignment Results")
    alignment_files = glob.glob("output/alignment/*.sam")
    
    for alignment_file in alignment_files:
        st.markdown(f"**{os.path.basename(alignment_file)}**")
        st.progress(100)  # Example progress bar (replace with actual progress logic if available)
        st.download_button(
            label="Download",
            data=open(alignment_file, "rb"),
            file_name=os.path.basename(alignment_file),
            mime="text/plain"
        )

# Tab 4: Variant Calling Results
with tabs[3]:
    st.subheader("Variant Calling Results")
    vcf_files = glob.glob("output/variants/*.vcf")
    
    for vcf_file in vcf_files:
        st.markdown(f"**{os.path.basename(vcf_file)}**")
        st.progress(100)  # Example progress bar (replace with actual progress logic if available)
        st.download_button(
            label="Download",
            data=open(vcf_file, "rb"),
            file_name=os.path.basename(vcf_file),
            mime="text/plain"
        )

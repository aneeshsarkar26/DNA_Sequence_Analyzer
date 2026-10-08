import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)


# ============================================================
# TITLE
# ============================================================

st.title("🧬 DNA Sequence Analyzer")

st.write(
    """
    A beginner-friendly bioinformatics application for analyzing
    DNA sequences, detecting mutations, searching motifs,
    translating DNA into protein, and visualizing sequence data.
    """
)


# ============================================================
# CODON TABLE
# ============================================================

CODON_TABLE = {

    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L",

    "UCU": "S", "UCC": "S",
    "UCA": "S", "UCG": "S",

    "UAU": "Y", "UAC": "Y",

    "UAA": "*",
    "UAG": "*",

    "UGU": "C", "UGC": "C",
    "UGA": "*",
    "UGG": "W",

    "CUU": "L", "CUC": "L",
    "CUA": "L", "CUG": "L",

    "CCU": "P", "CCC": "P",
    "CCA": "P", "CCG": "P",

    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",

    "CGU": "R", "CGC": "R",
    "CGA": "R", "CGG": "R",

    "AUU": "I", "AUC": "I",
    "AUA": "I",

    "AUG": "M",

    "ACU": "T", "ACC": "T",
    "ACA": "T", "ACG": "T",

    "AAU": "N", "AAC": "N",

    "AAA": "K", "AAG": "K",

    "AGU": "S", "AGC": "S",

    "AGA": "R", "AGG": "R",

    "GUU": "V", "GUC": "V",
    "GUA": "V", "GUG": "V",

    "GCU": "A", "GCC": "A",
    "GCA": "A", "GCG": "A",

    "GAU": "D", "GAC": "D",

    "GAA": "E", "GAG": "E",

    "GGU": "G", "GGC": "G",
    "GGA": "G", "GGG": "G"
}


# ============================================================
# AMINO ACID NAMES
# ============================================================

AMINO_ACIDS = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "E": "Glutamic acid",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine"
}


# ============================================================
# DNA COMPLEMENT
# ============================================================

def get_complement(sequence):

    translation_table = str.maketrans(
        "ATGC",
        "TACG"
    )

    return sequence.translate(translation_table)


# ============================================================
# REVERSE COMPLEMENT
# ============================================================

def get_reverse_complement(sequence):

    return get_complement(sequence)[::-1]


# ============================================================
# DNA → RNA
# ============================================================

def dna_to_rna(sequence):

    return sequence.replace("T", "U")


# ============================================================
# RNA → PROTEIN
# ============================================================

def rna_to_protein(rna, frame=0):

    protein = ""

    for i in range(
        frame,
        len(rna) - 2,
        3
    ):

        codon = rna[i:i + 3]

        amino_acid = CODON_TABLE.get(codon)

        if amino_acid is None:
            continue

        # Stop codon
        if amino_acid == "*":
            break

        protein += amino_acid

    return protein


# ============================================================
# FIND START CODONS
# ============================================================

def find_start_codons(sequence):

    positions = []

    start = 0

    while True:

        position = sequence.find(
            "ATG",
            start
        )

        if position == -1:
            break

        positions.append(position + 1)

        start = position + 1

    return positions


# ============================================================
# FIND STOP CODONS
# ============================================================

def find_stop_codons(sequence):

    stop_codons = ["TAA", "TAG", "TGA"]

    results = []

    for codon in stop_codons:

        start = 0

        while True:

            position = sequence.find(
                codon,
                start
            )

            if position == -1:
                break

            results.append(
                (
                    codon,
                    position + 1
                )
            )

            start = position + 1

    return sorted(
        results,
        key=lambda x: x[1]
    )


# ============================================================
# FIND MOTIF
# ============================================================

def find_motif(sequence, motif):

    positions = []

    start = 0

    while True:

        position = sequence.find(
            motif,
            start
        )

        if position == -1:
            break

        positions.append(position + 1)

        start = position + 1

    return positions


# ============================================================
# FIND ORFs
# ============================================================

def find_orfs(sequence):

    orfs = []

    stop_codons = {
        "TAA",
        "TAG",
        "TGA"
    }

    for frame in range(3):

        i = frame

        while i < len(sequence) - 2:

            codon = sequence[i:i + 3]

            if codon == "ATG":

                j = i + 3

                while j < len(sequence) - 2:

                    stop = sequence[j:j + 3]

                    if stop in stop_codons:

                        orf_sequence = sequence[
                            i:j + 3
                        ]

                        orfs.append({
                            "Frame": frame + 1,
                            "Start": i + 1,
                            "End": j + 3,
                            "Length": len(orf_sequence),
                            "Sequence": orf_sequence
                        })

                        break

                    j += 3

            i += 3

    return orfs


# ============================================================
# CLEAN DNA SEQUENCE
# ============================================================

def clean_sequence(sequence):

    sequence = sequence.upper()

    sequence = sequence.replace(
        " ",
        ""
    )

    sequence = sequence.replace(
        "\n",
        ""
    )

    sequence = sequence.replace(
        "\r",
        ""
    )

    return sequence


# ============================================================
# VALIDATE DNA
# ============================================================

def validate_dna(sequence):

    if sequence == "":
        return False

    for base in sequence:

        if base not in "ATGC":
            return False

    return True


# ============================================================
# GC CONTENT
# ============================================================

def calculate_gc(sequence):

    if len(sequence) == 0:
        return 0

    return (
        (
            sequence.count("G")
            +
            sequence.count("C")
        )
        /
        len(sequence)
    ) * 100


# ============================================================
# BASE COUNTS
# ============================================================

def get_base_counts(sequence):

    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }


# ============================================================
# SEQUENCE SIMILARITY
# ============================================================

def calculate_similarity(seq1, seq2):

    minimum_length = min(
        len(seq1),
        len(seq2)
    )

    if minimum_length == 0:
        return 0

    matches = 0

    for i in range(minimum_length):

        if seq1[i] == seq2[i]:

            matches += 1

    return (
        matches /
        minimum_length
    ) * 100


# ============================================================
# MUTATION DETECTION
# ============================================================

def detect_mutations(
    normal,
    mutated
):

    mutations = []

    max_length = max(
        len(normal),
        len(mutated)
    )

    for i in range(max_length):

        normal_base = (
            normal[i]
            if i < len(normal)
            else "-"
        )

        mutated_base = (
            mutated[i]
            if i < len(mutated)
            else "-"
        )

        if normal_base != mutated_base:

            if normal_base == "-":

                mutation_type = "Insertion"

            elif mutated_base == "-":

                mutation_type = "Deletion"

            else:

                mutation_type = "Substitution"

            mutations.append({
                "Position": i + 1,
                "Normal": normal_base,
                "Mutated": mutated_base,
                "Type": mutation_type
            })

    return mutations


# ============================================================
# FASTA PARSER
# ============================================================

def parse_fasta(text):

    sequences = {}

    current_name = None
    current_sequence = []

    for line in text.splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if current_name is not None:

                sequences[current_name] = (
                    "".join(current_sequence)
                )

            current_name = line[1:].strip()

            current_sequence = []

        else:

            current_sequence.append(
                line
            )

    if current_name is not None:

        sequences[current_name] = (
            "".join(current_sequence)
        )

    return sequences


# ============================================================
# SIDEBAR
# ============================================================




# ============================================================
# TABS
# ============================================================

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🧬 DNA Analyzer",
    "🧪 Protein Translation",
    "🔬 Mutation Analysis",
    "🔍 Motif & ORF",
    "📁 FASTA & Comparison"
])


# ============================================================
# TAB 1 — DNA ANALYZER
# ============================================================

with tab1:

    st.header(
        "🧬 DNA Sequence Analysis"
    )

    sequence_input = st.text_area(
        "Enter DNA sequence:",
        placeholder="Example: ATGCGTAGCTAGCTAGCGATCG",
        height=150
    )

    analyze_button = st.button(
        "🔬 Analyze DNA",
        type="primary"
    )


    if analyze_button:

        sequence = clean_sequence(
            sequence_input
        )

        if not sequence:

            st.warning(
                "⚠️ Please enter a DNA sequence."
            )

        elif not validate_dna(sequence):

            st.error(
                "❌ Invalid DNA sequence."
            )

            st.write(
                "Only A, T, G and C are allowed."
            )

        else:

            st.success(
                "✅ Valid DNA sequence!"
            )


            # ------------------------------------------------
            # BASIC ANALYSIS
            # ------------------------------------------------

            length = len(sequence)

            counts = get_base_counts(
                sequence
            )

            gc = calculate_gc(
                sequence
            )

            at = 100 - gc


            # ------------------------------------------------
            # METRICS
            # ------------------------------------------------

            col1, col2, col3 = st.columns(3)

            with col1:

                st.metric(
                    "Sequence Length",
                    f"{length} bp"
                )

            with col2:

                st.metric(
                    "GC Content",
                    f"{gc:.2f}%"
                )

            with col3:

                st.metric(
                    "AT Content",
                    f"{at:.2f}%"
                )


            # ------------------------------------------------
            # BASE COUNTS
            # ------------------------------------------------

            st.subheader(
                "🔢 Base Counts"
            )

            col1, col2, col3, col4 = st.columns(4)

            with col1:
                st.metric(
                    "A",
                    counts["A"]
                )

            with col2:
                st.metric(
                    "T",
                    counts["T"]
                )

            with col3:
                st.metric(
                    "G",
                    counts["G"]
                )

            with col4:
                st.metric(
                    "C",
                    counts["C"]
                )


            # ------------------------------------------------
            # PERCENTAGES
            # ------------------------------------------------

            st.subheader(
                "📊 Base Percentages"
            )

            percentage_data = {

                "Base": [
                    "A",
                    "T",
                    "G",
                    "C"
                ],

                "Count": [
                    counts["A"],
                    counts["T"],
                    counts["G"],
                    counts["C"]
                ]
            }

            df = pd.DataFrame(
                percentage_data
            )

            df["Percentage"] = (
                df["Count"] /
                length
            ) * 100

            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True
            )


            # ------------------------------------------------
            # SEQUENCE
            # ------------------------------------------------

            st.subheader(
                "🧬 Original Sequence"
            )

            st.code(
                sequence
            )


            # ------------------------------------------------
            # COMPLEMENT
            # ------------------------------------------------

            complement = get_complement(
                sequence
            )

            st.subheader(
                "🔄 Complement"
            )

            st.code(
                complement
            )


            # ------------------------------------------------
            # REVERSE COMPLEMENT
            # ------------------------------------------------

            reverse_complement = (
                get_reverse_complement(
                    sequence
                )
            )

            st.subheader(
                "🔁 Reverse Complement"
            )

            st.code(
                reverse_complement
            )


            # ------------------------------------------------
            # GRAPH
            # ------------------------------------------------

            st.subheader(
                "📈 Base Composition"
            )

            fig, ax = plt.subplots()

            ax.bar(
                ["A", "T", "G", "C"],
                [
                    counts["A"],
                    counts["T"],
                    counts["G"],
                    counts["C"]
                ]
            )

            ax.set_xlabel(
                "DNA Base"
            )

            ax.set_ylabel(
                "Count"
            )

            ax.set_title(
                "DNA Base Composition"
            )

            st.pyplot(fig)


# ============================================================
# TAB 2 — PROTEIN TRANSLATION
# ============================================================

with tab2:

    st.header(
        "🧪 DNA → RNA → Protein"
    )

    sequence_input = st.text_area(
        "Enter DNA sequence:",
        placeholder="Example: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
        key="protein_input"
    )


    if st.button(
        "🧬 Translate DNA",
        type="primary"
    ):

        sequence = clean_sequence(
            sequence_input
        )

        if not validate_dna(sequence):

            st.error(
                "❌ Invalid DNA sequence."
            )

        else:

            # DNA → RNA

            rna = dna_to_rna(
                sequence
            )

            st.subheader(
                "DNA"
            )

            st.code(
                sequence
            )


            st.subheader(
                "RNA"
            )

            st.code(
                rna
            )


            # --------------------------------------------
            # THREE READING FRAMES
            # --------------------------------------------

            st.subheader(
                "🧬 Reading Frames"
            )

            for frame in range(3):

                protein = rna_to_protein(
                    rna,
                    frame
                )

                st.write(
                    f"Reading Frame {frame + 1}"
                )

                if protein:

                    st.code(
                        protein
                    )

                else:

                    st.info(
                        "No protein could be translated."
                    )


            # --------------------------------------------
            # START CODONS
            # --------------------------------------------

            starts = find_start_codons(
                sequence
            )

            st.subheader(
                "▶️ Start Codons"
            )

            if starts:

                st.write(
                    "ATG found at positions:",
                    starts
                )

            else:

                st.write(
                    "No ATG start codon found."
                )


            # --------------------------------------------
            # STOP CODONS
            # --------------------------------------------

            stops = find_stop_codons(
                sequence
            )

            st.subheader(
                "⏹️ Stop Codons"
            )

            if stops:

                for codon, position in stops:

                    st.write(
                        f"{codon} → position {position}"
                    )

            else:

                st.write(
                    "No stop codons found."
                )


# ============================================================
# TAB 3 — MUTATION ANALYSIS
# ============================================================

with tab3:

    st.header(
        "🔬 DNA Mutation Detector"
    )

    st.write(
        """
        Enter a normal/reference DNA sequence and a mutated
        DNA sequence to compare them.
        """
    )


    normal_input = st.text_area(
        "Normal DNA sequence:",
        key="normal_sequence"
    )


    mutated_input = st.text_area(
        "Mutated DNA sequence:",
        key="mutated_sequence"
    )


    if st.button(
        "🔬 Detect Mutations",
        type="primary"
    ):

        normal = clean_sequence(
            normal_input
        )

        mutated = clean_sequence(
            mutated_input
        )


        if not validate_dna(normal):

            st.error(
                "❌ Invalid normal DNA sequence."
            )


        elif not validate_dna(mutated):

            st.error(
                "❌ Invalid mutated DNA sequence."
            )


        else:

            mutations = detect_mutations(
                normal,
                mutated
            )


            similarity = calculate_similarity(
                normal,
                mutated
            )


            st.metric(
                "Sequence Similarity",
                f"{similarity:.2f}%"
            )


            if mutations:

                st.subheader(
                    "🧬 Detected Mutations"
                )

                mutation_df = pd.DataFrame(
                    mutations
                )

                st.dataframe(
                    mutation_df,
                    use_container_width=True,
                    hide_index=True
                )


                st.write(
                    f"Total differences: "
                    f"**{len(mutations)}**"
                )

            else:

                st.success(
                    "✅ No differences detected."
                )


# ============================================================
# TAB 4 — MOTIF & ORF ANALYSIS
# ============================================================

with tab4:

    st.header(
        "🔍 Motif & ORF Analysis"
    )


    sequence_input = st.text_area(
        "Enter DNA sequence:",
        key="motif_sequence"
    )


    motif_input = st.text_input(
        "Enter motif to search:",
        placeholder="Example: ATG",
        key="motif"
    )


    if st.button(
        "🔍 Search",
        type="primary"
    ):

        sequence = clean_sequence(
            sequence_input
        )

        motif = clean_sequence(
            motif_input
        )


        if not validate_dna(sequence):

            st.error(
                "❌ Invalid DNA sequence."
            )


        elif not motif:

            st.warning(
                "Please enter a motif."
            )


        elif not validate_dna(motif):

            st.error(
                "❌ Invalid motif."
            )


        else:

            # --------------------------------------------
            # MOTIF SEARCH
            # --------------------------------------------

            positions = find_motif(
                sequence,
                motif
            )


            st.subheader(
                "🔍 Motif Search"
            )


            if positions:

                st.success(
                    f"Motif found {len(positions)} time(s)."
                )

                st.write(
                    "Positions:",
                    positions
                )

            else:

                st.info(
                    "Motif not found."
                )


            # --------------------------------------------
            # ORF SEARCH
            # --------------------------------------------

            st.subheader(
                "🧬 Open Reading Frames"
            )


            orfs = find_orfs(
                sequence
            )


            if orfs:

                orf_df = pd.DataFrame(
                    orfs
                )

                st.dataframe(
                    orf_df,
                    use_container_width=True,
                    hide_index=True
                )

            else:

                st.info(
                    "No complete ORFs were found."
                )


# ============================================================
# TAB 5 — FASTA & SEQUENCE COMPARISON
# ============================================================

with tab5:

    st.header(
        "📁 FASTA File & Sequence Comparison"
    )


    # ========================================================
    # FASTA UPLOAD
    # ========================================================

    st.subheader(
        "📁 Upload FASTA File"
    )


    uploaded_file = st.file_uploader(
        "Upload a .fasta or .fa file:",
        type=["fasta", "fa", "txt"]
    )


    if uploaded_file is not None:

        fasta_text = (
            uploaded_file
            .getvalue()
            .decode("utf-8")
        )


        sequences = parse_fasta(
            fasta_text
        )


        if sequences:

            st.success(
                f"Loaded {len(sequences)} sequence(s)."
            )


            for name, sequence in sequences.items():

                sequence = clean_sequence(
                    sequence
                )


                st.write(
                    f"**{name}**"
                )

                st.code(
                    sequence
                )


                if validate_dna(sequence):

                    counts = get_base_counts(
                        sequence
                    )

                    gc = calculate_gc(
                        sequence
                    )


                    col1, col2, col3 = st.columns(3)


                    with col1:

                        st.metric(
                            "Length",
                            f"{len(sequence)} bp"
                        )


                    with col2:

                        st.metric(
                            "GC Content",
                            f"{gc:.2f}%"
                        )


                    with col3:

                        st.metric(
                            "A + T",
                            counts["A"] +
                            counts["T"]
                        )


                else:

                    st.warning(
                        "This sequence contains invalid DNA bases."
                    )


        else:

            st.error(
                "Could not read the FASTA file."
            )


    # ========================================================
    # TWO-SEQUENCE COMPARISON
    # ========================================================

    st.divider()

    st.subheader(
        "🔬 Compare Two DNA Sequences"
    )


    seq1_input = st.text_area(
        "Sequence 1:",
        key="comparison_seq1"
    )


    seq2_input = st.text_area(
        "Sequence 2:",
        key="comparison_seq2"
    )


    if st.button(
        "📊 Compare Sequences"
    ):

        seq1 = clean_sequence(
            seq1_input
        )

        seq2 = clean_sequence(
            seq2_input
        )


        if not validate_dna(seq1):

            st.error(
                "❌ Sequence 1 is invalid."
            )


        elif not validate_dna(seq2):

            st.error(
                "❌ Sequence 2 is invalid."
            )


        else:

            similarity = calculate_similarity(
                seq1,
                seq2
            )


            st.metric(
                "Sequence Similarity",
                f"{similarity:.2f}%"
            )


            st.write(
                f"Sequence 1 length: "
                f"**{len(seq1)} bp**"
            )

            st.write(
                f"Sequence 2 length: "
                f"**{len(seq2)} bp**"
            )


# ============================================================
# FOOTER
# ============================================================

st.divider()

st.caption(
    "🧬 DNA Sequence Analyzer | "
    "Educational Bioinformatics Project"
)

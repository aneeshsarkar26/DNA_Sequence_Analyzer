import streamlit as st
import matplotlib.pyplot as plt


# --------------------------------------------------
# PAGE CONFIGURATION
# --------------------------------------------------

st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)


# --------------------------------------------------
# TITLE
# --------------------------------------------------

st.title("🧬 DNA Sequence Analyzer")

st.write(
    "Enter a DNA sequence below and the program will analyze "
    "its composition and convert it into RNA and protein."
)


# --------------------------------------------------
# CODON TABLE
# --------------------------------------------------

codon_table = {
    # Phenylalanine
    "UUU": "F",
    "UUC": "F",

    # Leucine
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",

    # Isoleucine
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",

    # Methionine / Start
    "AUG": "M",

    # Valine
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",

    # Serine
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "AGU": "S",
    "AGC": "S",

    # Proline
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",

    # Threonine
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",

    # Alanine
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",

    # Tyrosine
    "UAU": "Y",
    "UAC": "Y",

    # Histidine
    "CAU": "H",
    "CAC": "H",

    # Glutamine
    "CAA": "Q",
    "CAG": "Q",

    # Asparagine
    "AAU": "N",
    "AAC": "N",

    # Lysine
    "AAA": "K",
    "AAG": "K",

    # Aspartic acid
    "GAU": "D",
    "GAC": "D",

    # Glutamic acid
    "GAA": "E",
    "GAG": "E",

    # Cysteine
    "UGU": "C",
    "UGC": "C",

    # Tryptophan
    "UGG": "W",

    # Arginine
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",

    # Glycine
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",

    # Stop codons
    "UAA": "*",
    "UAG": "*",
    "UGA": "*"
}


# --------------------------------------------------
# DNA COMPLEMENT FUNCTION
# --------------------------------------------------

def get_complement(sequence):

    complement = ""

    for base in sequence:

        if base == "A":
            complement += "T"

        elif base == "T":
            complement += "A"

        elif base == "G":
            complement += "C"

        elif base == "C":
            complement += "G"

    return complement


# --------------------------------------------------
# PROTEIN TRANSLATION FUNCTION
# --------------------------------------------------

def translate_rna(rna):

    protein = ""

    # Read RNA in groups of 3
    for i in range(0, len(rna) - 2, 3):

        codon = rna[i:i + 3]

        # Check whether codon exists
        if codon in codon_table:

            amino_acid = codon_table[codon]

            # Stop translation at stop codon
            if amino_acid == "*":
                break

            protein += amino_acid

    return protein


# --------------------------------------------------
# DNA INPUT
# --------------------------------------------------

sequence_input = st.text_area(
    "Enter your DNA sequence:",
    placeholder="Example: ATGCGTAGCTAGCTAGCGATCG",
    height=120
)


# --------------------------------------------------
# ANALYZE BUTTON
# --------------------------------------------------

if st.button("🔬 ANALYZE DNA", type="primary"):

    # Convert to uppercase
    sequence = sequence_input.upper()

    # Remove spaces and new lines
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    sequence = sequence.replace("\r", "")

    # --------------------------------------------------
    # CHECK IF USER ENTERED SOMETHING
    # --------------------------------------------------

    if sequence == "":

        st.warning("⚠️ Please enter a DNA sequence.")

    else:

        # --------------------------------------------------
        # VALIDATE DNA
        # --------------------------------------------------

        valid = True

        for base in sequence:

            if base not in "ATGC":

                valid = False
                break


        # --------------------------------------------------
        # INVALID DNA
        # --------------------------------------------------

        if not valid:

            st.error(
                "❌ Invalid DNA sequence! "
                "Only A, T, G and C are allowed."
            )


        # --------------------------------------------------
        # VALID DNA
        # --------------------------------------------------

        else:

            st.success("✅ Valid DNA sequence!")


            # --------------------------------------------------
            # BASIC CALCULATIONS
            # --------------------------------------------------

            length = len(sequence)

            a_count = sequence.count("A")
            t_count = sequence.count("T")
            g_count = sequence.count("G")
            c_count = sequence.count("C")


            # --------------------------------------------------
            # GC AND AT CONTENT
            # --------------------------------------------------

            gc_content = ((g_count + c_count) / length) * 100

            at_content = ((a_count + t_count) / length) * 100


            # --------------------------------------------------
            # BASIC INFORMATION
            # --------------------------------------------------

            st.header("📊 Sequence Analysis")


            # First row of metrics

            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric(
                    "Sequence Length",
                    f"{length} bp"
                )

            with col2:
                st.metric(
                    "GC Content",
                    f"{gc_content:.2f}%"
                )

            with col3:
                st.metric(
                    "AT Content",
                    f"{at_content:.2f}%"
                )


            # --------------------------------------------------
            # BASE COUNTS
            # --------------------------------------------------

            st.subheader("🧬 Base Composition")

            col1, col2, col3, col4 = st.columns(4)

            with col1:
                st.metric("Adenine (A)", a_count)

            with col2:
                st.metric("Thymine (T)", t_count)

            with col3:
                st.metric("Guanine (G)", g_count)

            with col4:
                st.metric("Cytosine (C)", c_count)


            # --------------------------------------------------
            # ORIGINAL SEQUENCE
            # --------------------------------------------------

            st.subheader("🧬 Original DNA Sequence")

            st.code(sequence)


            # --------------------------------------------------
            # COMPLEMENT
            # --------------------------------------------------

            complement = get_complement(sequence)

            st.subheader("🔄 Complementary DNA")

            st.code(complement)


            # --------------------------------------------------
            # REVERSE COMPLEMENT
            # --------------------------------------------------

            reverse_complement = complement[::-1]

            st.subheader("🔁 Reverse Complement")

            st.code(reverse_complement)


            # --------------------------------------------------
            # DNA → RNA
            # --------------------------------------------------

            rna = sequence.replace("T", "U")

            st.subheader("🧪 DNA → RNA")

            st.code(rna)


            # --------------------------------------------------
            # RNA → PROTEIN
            # --------------------------------------------------

            protein = translate_rna(rna)

            st.subheader("🧬 RNA → Protein")

            if protein:

                st.code(protein)

            else:

                st.info(
                    "No complete codon could be translated."
                )


            # --------------------------------------------------
            # GRAPH
            # --------------------------------------------------

            st.header("📈 Base Composition Graph")

            bases = ["A", "T", "G", "C"]

            counts = [
                a_count,
                t_count,
                g_count,
                c_count
            ]


            # Create graph

            fig, ax = plt.subplots()

            ax.bar(
                bases,
                counts
            )

            ax.set_xlabel("DNA Base")

            ax.set_ylabel("Number of Bases")

            ax.set_title(
                "DNA Base Composition"
            )


            # Display graph

            st.pyplot(fig)


            # --------------------------------------------------
            # SUMMARY
            # --------------------------------------------------

            st.header("📋 Analysis Summary")

            st.write(
                f"""
                **Sequence Length:** {length} base pairs

                **Adenine (A):** {a_count}

                **Thymine (T):** {t_count}

                **Guanine (G):** {g_count}

                **Cytosine (C):** {c_count}

                **GC Content:** {gc_content:.2f}%

                **AT Content:** {at_content:.2f}%
                """
            )
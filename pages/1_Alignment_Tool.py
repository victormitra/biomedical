import streamlit as st
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Alignment metrics
def compute_metrics(seq1, seq2):
    aligned_length = len(seq1)
    match_count = sum(a == b for a, b in zip(seq1, seq2))
    percent_identity = (match_count / aligned_length) * 100 if aligned_length > 0 else 0
    return match_count, aligned_length, percent_identity

# Pretty text alignment
def pretty_alignment(seq1, seq2, width=60):
    lines = []
    for i in range(0, len(seq1), width):
        s1 = seq1[i:i+width]
        s2 = seq2[i:i+width]
        match_line = ''.join('|' if a == b else ' ' for a, b in zip(s1, s2))

        lines.append(f"Seq1: {s1}")
        lines.append(f"      {match_line}")
        lines.append(f"Seq2: {s2}")
        lines.append("")
    return "\n".join(lines)

# Streamlit page UI
st.title("🧬 Protein Sequence Aligner")

seq1 = st.text_area("Enter Sequence 1", value="MYGKIIFVLLLSAIVSISASSTTGV")
seq2 = st.text_area("Enter Sequence 2", value="MYGKIIFVLAASTTGVAMHTST")

method = st.selectbox("Alignment Method", ["Needleman-Wunsch (Global)", "Smith-Waterman (Local)"])

if st.button("Align"):
    if "Global" in method:
        alignments = pairwise2.align.globalxx(seq1, seq2)
    else:
        alignments = pairwise2.align.localxx(seq1, seq2)

    align1, align2, score, start, end = alignments[0]

    # Alignment Display
    st.subheader("Alignment")
    st.code(pretty_alignment(align1, align2), language="text")

    # Metrics
    matches, aligned_len, identity = compute_metrics(align1, align2)
    st.subheader("Metrics")
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Score", f"{score:.2f}")
    col2.metric("Matches", matches)
    col3.metric("Aligned Length", aligned_len)
    col4.metric("Percent Identity", f"{identity:.2f}%")

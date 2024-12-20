function highlightNucleotides() {
    // Get the DNA sequence element
    const sequenceElement = document.getElementById("dna-sequence");
    
    // Get plain text content (ignoring existing tags)
    let sequence = sequenceElement.textContent;

    // Replace each nucleotide with styled spans
    sequence = sequence
        .replace(/A/g, "<span class='A'>A</span>")
        .replace(/T/g, "<span class='T'>T</span>")
        .replace(/C/g, "<span class='C'>C</span>")
        .replace(/G/g, "<span class='G'>G</span>");

    // Update the element's innerHTML
    sequenceElement.innerHTML = sequence;
}

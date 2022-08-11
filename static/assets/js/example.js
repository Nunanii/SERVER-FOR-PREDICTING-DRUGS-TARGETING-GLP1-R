var content = "SMILES\nCC1=C(N2CCCC(=O)NC2=C1C#N)C\nc1ccccc1\nc\nCC(=O)Nc1ccc(cc1)O\nC1=CC=C(C=C1)CON=CC2=CC=CC=N2";
const example_file = document.getElementById("example_file");
const blob1 = new Blob([content], {type: "text/plain"})
example_file.href = URL.createObjectURL(blob1)


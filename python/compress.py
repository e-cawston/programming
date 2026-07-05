import fitz  # PyMuPDF

def compress_pdf_no_image_changes(input_pdf, output_pdf):
    doc = fitz.open(input_pdf)

    # Save without image recompression
    doc.save(
        output_pdf,
        deflate=True,   # compress text/objects only
        clean=True,     # remove unused objects
        garbage=1       # light garbage collection (prevents freezing)
    )

    doc.close()

# Put your filenames here:
compress_pdf_no_image_changes("Euan_Cawston_status-funding.pdf", "Euan_Cawston_status-funding_compressed.pdf")

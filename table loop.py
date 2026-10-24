from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from prettytable import PrettyTable


def filter_hits(blast_record, evalue_threshold=1e-5, query_cover_range=(40, 100), identity_range=(70, 100), score_threshold=200):
    filtered_hits = []

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            
            if (
                hsp.expect <= evalue_threshold and
                query_cover_range[0] <= (hsp.query_end - hsp.query_start + 1) / len(query_sequence) * 100 <= query_cover_range[1] and
                identity_range[0] <= (hsp.identities / hsp.align_length) * 100 <= identity_range[1] and
                hsp.score > score_threshold
            ):
                # Extract species name from the alignment title
                species_name = alignment.title.split('|')[-1].strip()
                
                percentage_identity = (hsp.identities / hsp.align_length) * 100
                query_cover = (hsp.query_end - hsp.query_start + 1) / len(query_sequence) * 100

                filtered_hits.append({
                    'Species Name': species_name,
                    'Accession': alignment.accession,
                    'Score': hsp.score,
                    'E-Value': hsp.expect,
                    'Query Cover': query_cover,
                    'Percentage Identity': percentage_identity
                })

    # Sort hits by E-value (ascending), query cover (descending), and percentage identity (descending)
    sorted_hits = sorted(filtered_hits, key=lambda x: (x['E-Value'], -x['Query Cover'], -x['Percentage Identity']))

    # Take the top three hits
    top_three_hits = sorted_hits[:3]

    return top_three_hits


#prog starts here
#variable declaration
query=[]
tmp=[]
msg=''

with open('input1.fasta','r') as file:
    for line in file:
        if line !='\n':
            tmp.append(line)
        if line =='\n':
            query.append(tmp)
            tmp=[]
    file.close()

for i in range(0,len(query)):
    query_sequence = query[i][1]
    result_handle = NCBIWWW.qblast("blastx", "refseq_protein", query_sequence)
    blast_record = NCBIXML.read(result_handle)
    top_three_hits = filter_hits(blast_record, evalue_threshold=1e-5, query_cover_range=(40, 100), identity_range=(70, 100), score_threshold=200)
    if not top_three_hits:
        print("No matching result")
        msg=str(query[i][0])
        table_file = open("top_three_blast_hits_table.txt", "a")
        table_file.write(str(msg))
        table_file.write(str('Has no result\n\n'))
    
    else:
        table = PrettyTable()
        table.field_names = ['Species Name','Accession', 'Score', 'E-Value', 'Query Cover', 'Percentage Identity']
        for hit in top_three_hits:
            table.add_row([hit['Species Name'],hit['Accession'], hit['Score'], hit['E-Value'], f"{hit['Query Cover']:.2f}%", f"{hit['Percentage Identity']:.2f}%"])
        print(table)
        with open("top_three_blast_hits_table.txt", "a") as table_file:
            table_file.write(str(table))
            table_file.write(str('\n'))
            table_file.write(str(query[i][0]))
            table_file.write(str('\n'))        
    result_handle.close()

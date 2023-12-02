from flask import Flask, request, jsonify
import random
from Bio import SeqIO

app = Flask(__name__)

def parse_fasta(fasta_content):
    reads = []
    current_read = None

    # Split the content into lines
    lines = fasta_content.split('\n')

    for line in lines:
        line = line.strip()

        # Check for a FASTA header line
        if line.startswith('>'):
            # Start a new read
            current_read = {'header': line[1:], 'sequence': ''}
            reads.append(current_read)
        elif current_read is not None:
            # Add sequence lines to the current read
            current_read['sequence'] += line

    return reads


def debruijnize(reads, ind, out):
      nodes = set()
      not_starts = set()
      edges = []
      for r in reads:
          r1 = r[:-1]
          r2 = r[1:]
          nodes.add(r1)
          nodes.add(r2)
          edges.append((r1,r2))
          if r1 in out:
            out[r1] = out[r1]+1
          else :
            out[r1] = 1

          if r2 in ind:
            ind[r2] = ind[r2]+1
          else :
            ind[r2] = 1
          not_starts.add(r2)
      return (nodes,edges,list(nodes-not_starts))


def make_graph(edges):
    node_edge_map = {}
    for e in edges:
        n = e[0]
        if n in node_edge_map:
            node_edge_map[n].append(e[1])
        else:
            node_edge_map[n] = [e[1]]
    return node_edge_map

def getResults(request):
  fasta_file = request.files['fasta_file']
  kmer_length = int(request.form['kmer_length'])
  fasta_content = fasta_file.read().decode('utf-8')
  sequences = parse_fasta(fasta_content)
  reads = []

  for sequence in sequences:
    reads.append(sequence['sequence'])
  
  k_mers = []
  for i in reads:
    for j in range(0,len(i)-kmer_length+1):
      temp = i[j:j+kmer_length]
      if len(temp)==kmer_length:
        k_mers.append(temp)


  ind = {}
  out = {}
  
  G = debruijnize(k_mers, ind, out)
  c = 0
  s = 0
  o = 0
  start_nodes = []
  for i in G[0]:
    if i not in ind:
      ind[i] = 0
    if i not in out:
      out[i] = 0

    if ind[i] == out[i]:
      start_nodes.append(i)
      c+=1
    if abs(ind[i]-out[i])>1:
      o+=1
    else:
      s+=1
      if out[i]-ind[i] == 1:
        start_nodes.append(i)

  if len(start_nodes)>100:
    start_nodes = random.sample(start_nodes,100)
  
  m = make_graph(G[1])
  results = []
  for start in start_nodes:
    path = []
    visited = {}
    for i in G[0]:
      visited[i] = False
      if i not in m:
        m[i] = []
    def eulerian_dfs(m,v):
      visited[v] = True
      for next in m[v]:
        if visited[next]==False:
          eulerian_dfs(m,next)

      path.append(v)

    eulerian_dfs(m,start)
    path.reverse()

    def assemble_trail(trail):
        if len(trail) == 0:
            return ""
        result = trail[0][:-1]
        for node in trail:
            result += node[-1]
        return result
    result_genome = assemble_trail(path)
    results.append(result_genome)

  return results

def global_alignment(seq1, seq2, match_score=1, mismatch_penalty=0, gap_penalty=0):
    # Initialize the score matrix
    rows, cols = len(seq1) + 1, len(seq2) + 1
    score_matrix = [[0] * cols for _ in range(rows)]


    # Initialize the traceback matrix
    traceback_matrix = [[0] * cols for _ in range(rows)]

    # Initialize the first row and column of the score matrix
    for i in range(1, rows):
        score_matrix[i][0] = i * gap_penalty
        traceback_matrix[i][0] = 1  # 1 represents a gap in the traceback matrix

    for j in range(1, cols):
        score_matrix[0][j] = j * gap_penalty
        traceback_matrix[0][j] = 2  # 2 represents a gap in the traceback matrix

    # Fill in the score matrix and traceback matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty

            # Choose the maximum score
            max_score = max(match, delete, insert)
            score_matrix[i][j] = max_score

            # Update the traceback matrix
            if max_score == match:
                traceback_matrix[i][j] = 0  # 0 represents a match in the traceback matrix
            elif max_score == delete:
                traceback_matrix[i][j] = 1  # 1 represents a gap in the traceback matrix
            else:
                traceback_matrix[i][j] = 2  # 2 represents a gap in the traceback matrix

    # Perform traceback to find the aligned sequences
    aligned_seq1, aligned_seq2 = "", ""
    i, j = rows - 1, cols - 1

    
    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 0:  # Match
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 1:  # Gap in seq1
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:  # Gap in seq2
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    return aligned_seq1,aligned_seq2,score_matrix[rows-1][cols-1],1,1


def getStats(results):
  mn = len(results[0])
  mx = len(results[0])
  sum = 0

  for result in results:
    length = len(result)
    mn = min(mn, length)
    mx = max(mx, length)
    sum = sum + length
  
  avg = sum / len(results)

  return mn, mx, avg
  

         


@app.route('/process_fasta', methods=['POST'])
def process_fasta():
  try:
    results = getResults(request)
    stats = getStats(results)
    return jsonify({'results': results, 'minLength': stats[0], 'maxLength': stats[1], 'avgLength': stats[2]}), 200
  except Exception as e:
    return jsonify({"error": str(e)}), 500
  


@app.route('/fasta_accuracy', methods=['POST'])
def fasta_accuracy():
  try:
    results = request.form.getlist('results')
    print('length of results array is: ',len(results))
    match_score = int(request.form['match_score'])
    print('match_score is: ', match_score)
    mismatch_penalty = int(request.form['mismatch_penalty'])
    print('mismatch_penalty is: ', mismatch_penalty)
    gap_penalty = int(request.form['gap_penalty'])
    print('gap_penalty is: ', gap_penalty)
    key_file = request.files['key_file']
    print('key file loaded successfully')
    genome = key_file.read().decode('utf-8')
    print('length of genome is: ' ,len(genome))

    allResults = []
    seq1 = ''
    seq2 = ''
    accuracy = 0
    best = results[0]
    mismatch = ''
    gap = ''

    for gene in results:
      if(1==1):
         allResults.append(1)
         continue
      answer = global_alignment(genome, gene, match_score, mismatch_penalty, gap_penalty)
      if(answer[2]>accuracy):
        accuracy = answer[2]
        seq1 = answer[0]
        seq2 = answer[1]
        best = gene
        mismatch = answer[3]
        gap = answer[4]
      
      temp = {
        'seq1': answer[0],
        'seq2': answer[1],
        'acc': answer[2]
      }

      allResults.append(temp)

      stats = getStats(results)
      
    return jsonify({'allResults': allResults, 'results': results, 'seq1': seq1, 'seq2': seq2, 'accuracy': accuracy,
                    'minLength': stats[0], 'maxLength': stats[1], 'avgLength': stats[2], 'best': best, 'mismatch': mismatch, 'gap': gap}), 200

  except Exception as e:
    return jsonify({"error": str(e)}), 500

  


if(__name__ == '__main__'):
  app.run()
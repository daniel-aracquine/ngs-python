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



@app.route('/process_fasta', methods=['POST'])
def process_fasta():
  try:
    print('sarthak')
    fasta_file = request.files['fasta_file']
    kmer_length = int(request.form['kmer_length'])
    print('samarth', + kmer_length)
    print(fasta_file)
    fasta_content = fasta_file.read().decode('utf-8')
    print(fasta_content)
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

    return jsonify({'results': results}), 200
  except Exception as e:
    return jsonify({"error": str(e)}), 500

if(__name__ == '__main__'):
  app.run()
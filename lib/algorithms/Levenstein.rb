#!/usr/bin/ruby19

require 'narray'

# Inspired by http://www.informit.com/articles/article.aspx?p=683059&seqNum=36
class String

  def levenshtein(other, ins=2, del=2, sub=1)
    # ins, del, sub are weighted costs
    return nil if self.nil?
    return nil if other.nil?
    dm = []        # distance matrix

    # Initialize first row values
    dm[0] = (0..self.length).collect { |i| i * ins }
    fill = [0] * (self.length - 1)

    # Initialize first column values
    for i in 1..other.length
      dm[i] = [i * del, fill.flatten]
    end

    # populate matrix
    for i in 1..other.length
      for j in 1..self.length
    # critical comparison
        dm[i][j] = [
             dm[i-1][j-1] +
               (self[j-1] == other[i-1] ? 0 : sub),
                 dm[i][j-1] + ins,
             dm[i-1][j] + del
       ].min
      end
    end

    # The last value in matrix is the
    # Levenshtein distance between the strings
    [ dm, dm[other.length][self.length] ]
  end

end

# New Similarity Feature: In the paper An Information-Theoretic Definition of Similarity (6), Dekang Lin, Department of Computer Science, University of Manitoba, has stated that edit distances can be converted into similarity values as follows:
#   simedit(x, y) = 1/(1 + editDist(x, y))
# http://www.miislita.com/searchito/levenshtein-edit-distance.html
# http://www.cs.ualberta.ca/~lindek/papers/sim.pdf
def simedit( x, y, edit_distance )
  1 / ( 1 + edit_distance( x, y ) )
end


def similar( edit_distance )
  similar = NArray.float( edit_distance.length, edit_distance.first.length ) 

  for i in 0..(edit_distance.length - 1) do
    for j in 0..(edit_distance.first.length - 1) do
      similar[ i, j ] = 1.0 / ( 1.0 + edit_distance[i][j] )
    end
  end

  # the bigger the number the more similar
  p similar
end


if __FILE__ == $0
  raise ArgumentError, "Needs two files to compare" unless( ARGV.length == 2 )

  file1 = ( File.open( ARGV[0].to_s, "r" ).readlines.collect { |l| (l.chomp)  } ).join( "" )
  file2 = ( File.open( ARGV[1].to_s, "r" ).readlines.collect { |l| (l.chomp ) } ).join( "" )

  distance_matrix, distance_value = file1.levenshtein( file2 )
  p distance_value
  
  # similar( distance_matrix )

end


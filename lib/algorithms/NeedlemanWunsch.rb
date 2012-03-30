#!/usr/bin/ruby19


# Necessary libraries
require 'narray'



# @class        class NeedlemanWunsch # {{{
# @brief        Needleman Wunsch Algorithm implementation in Ruby
class NeedlemanWunsch

  # @fn         def initialize # {{{
  # @brief      Constructor of NeedlemanWunsch class
  #
  # @param      [Logger]        logger    Instantiated Logger class
  # @param      [Integer]       gap       Gap penalty of when introducing new gaps into the sequences (constant)
  def initialize logger = nil, genetic_sequence = false

    @logger               = logger
    @gap, @similarity     = ( genetic_sequence ) ? ( init_genetic ) : ( init_generic )

  end


  # @fn         def init_genetic # {{{
  # @brief      Initializes class with similarity matrix and gap penalty suitable for genetic sequences (Bio-Informatics)
  #
  # FIXME:      These values should be mapped out to a config file etc.
  def init_genetic

    gap = -5

    # similarity matrix when using genetic sequences
    #
    #     A  G  C  T
    #  A 10 -1 -3 -4
    #  G -1  7 -5 -3
    #  C -3 -5  9  0
    #  T -4 -3  0  8

    similarity = {  'AA' => 10,
                    'AG' => -1,
                    'AC' => -3,
                    'AT' => -4,
                    'GA' => -1,
                    'GG' =>  7,
                    'GC' => -5,
                    'GT' => -3,
                    'CA' => -3,
                    'CG' => -5,
                    'CC' =>  9,
                    'CT' =>  0,
                    'TA' => -4,
                    'TG' => -3,
                    'TC' =>  0,
                    'TT' =>  8 }

    return [ gap, similarity ]
  end # of def init_genetic # }}}


  # @fn         def init_generic # {{{
  # @brief      Initializes class with similarity matrix and gap penalty suitable for generic sequences (found outside of Bio-Informatics)
  #
  # FIXME:      These values should be mapped out to a config file etc.
  def init_generic

    gap = -2

    similarity = {  'AA' => 10,
                    'AG' => -1,
                    'AC' => -3,
                    'AT' => -4,
                    'GA' => -1,
                    'GG' =>  7,
                    'GC' => -5,
                    'GT' => -3,
                    'CA' => -3,
                    'CG' => -5,
                    'CC' =>  9,
                    'CT' =>  0,
                    'TA' => -4,
                    'TG' => -3,
                    'TC' =>  0,
                    'TT' =>  8 }

    return [ gap, similarity ]
  end # of def init_genetic # }}}


  def needle sequence, reference

    gap = @gap
    d = @similarity

    s = {}

    # All permutations e.g. 01, 10 etc. are penalized
    ["0", "1", "2", "3", "4", "5", "6", "7"].permutation(2) do |x|
      s[ x.join( "" ).to_s ] = -1
    end

    # Matches "00" and same letter matches .. similar are awarded
    ["0", "1", "2", "3", "4", "5", "6", "7"].each do |x|
      s[ ( x + x ).to_s ] = 1
    end

    rows = reference.length + 1
    cols = sequence.length + 1

    a = NArray.int( rows, cols )

    for i in 0...(rows) do a[i,0] = 0 end
    for j in 0...(cols) do a[0,j] = 0 end
    for i in 1...(rows)
      for j in 1...(cols)
        sec = s[ (reference[i-1].chr + sequence[j-1].chr) ]
        sec = 0 if( sec.nil? )
        choice1 = a[i-1, j-1] + sec
        # choice1 = a[i-1, j-1] + s[ (reference[i-1].chr + sequence[j-1].chr) ]

        choice2 = a[i-1, j] + gap
        choice3 = a[i, j-1] + gap
        a[i,j] = [choice1, choice2, choice3].max
      end
    end

    ref = ''
    seq = ''

    i = reference.length 
    j = sequence.length

    while (i > 0 and j > 0)
      score = a[i,j]
      score_diag = a[i-1, j-1]
      score_up = a[i, j-1]
      score_left = a[i-1, j]

      sec = s[ (reference[i-1].chr + sequence[j-1].chr) ]
      sec = 0 if( sec.nil? )

      if (score == score_diag + sec )
        ref = reference[i-1].chr + ref
        seq = sequence[j-1].chr + seq
        i -= 1
        j -= 1
      elsif (score == score_left + gap)
        ref = reference[i-1].chr + ref
        seq = '-' + seq
        i -= 1
      elsif (score == score_up + gap)
        ref = '-' + ref
        seq = sequence[j-1].chr + seq
        j -= 1
      end
    end

    while (i > 0)
      ref = reference[i-1].chr + ref
      seq = '-' + seq
      i -= 1    
    end

    while (j > 0)
      ref = '-' + ref
      seq = sequence[j-1].chr + seq
      j -= 1
    end

    [seq, ref] 
  end


end # of class NeedlemanWunsch # }}}




if __FILE__ == $0

  nw          = NeedlemanWunsch.new( nil, true )

  sequence    = "gccctagcg"
  reference   = "gcgcaatg"

  s, r =  nw.needle(sequence, reference)
  ns   = []
  0.upto( s.length - 1 ).each do |n|
    ns << n.to_s
  end

  ss = ns.zip( s.split( "" ) )
  ss.collect!{ |array| array.join( " " ) }
  ss = ss.join( "\n" )

  rs = ns.zip( r.split( "" ) )
  rs.collect!{ |array| array.join( " " ) }
  rs = rs.join( "\n" )


  #puts ss
  puts rs

end



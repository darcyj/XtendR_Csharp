using System;
using System.IO;
using System.Linq;
using System.Diagnostics;
using System.Text.RegularExpressions;
using System.Collections.Generic;
using System.Threading.Tasks;
using System.Collections.Concurrent;

namespace XtendR_Csharp
{
	class Functions{

	// OutputIndenter
	// indents text for output and log purposes
		public string OutputIndenter(string text, int indent_lvl, char indent=' ', int multiplier=4){
			string pre = new String (indent, indent_lvl * multiplier);
			string output = pre + text;
			return output;
		}

	// FastaReader
	// reads fasta file into array of strings
		public string[] FastaReader (string filepath, bool noisy=true, int indent_lvl=0){
			if(noisy){
				// Console.WriteLine ("Reading sequence data from " + filepath + ":");
				string tp0 = "Reading sequence data from " + filepath + ":";
				Console.WriteLine (OutputIndenter (tp0, indent_lvl));
			}
			Stopwatch watch = new Stopwatch();
			watch.Start ();
			string[] lines = File.ReadAllLines(filepath).Where((s) => !s.StartsWith(">")).ToArray();
			watch.Stop ();
			if(noisy){
				string time = watch.Elapsed.ToString(@"m\:ss");
				string tp = "Read " + lines.Length + " sequences in " + time + " minutes.";
				Console.WriteLine (OutputIndenter (tp, indent_lvl + 1));
			}
			return lines;
		}

	//ParallelKmerSearch
	// searches inputseqs for kmer - returns array of matching seqs
		public string[] ParallelKmerSearch (string[] inputseqs, string kmer, int ncores = 4, bool noisy = true, int indent_lvl=0, string input_desc = "raw seqs"){
			// tell user what we're doing
			if (noisy) {
				string tp = "Searching " + input_desc + " for kmer " + kmer;
				Console.WriteLine (OutputIndenter (tp, indent_lvl));
			}
			Stopwatch watch = new Stopwatch();
			watch.Start ();

			// adding elements to a list is faster than concatenating to an array
			// But Lists aren't thread-safe I guess. So I have to use a ConcurrentQueue instead.
			ConcurrentQueue<string> matches = new ConcurrentQueue<string>();

			// regex object for kmer
			Regex r = new Regex(kmer, RegexOptions.None);

			// for each string, look for matches
			Parallel.ForEach (
				inputseqs, 
				new ParallelOptions { MaxDegreeOfParallelism = ncores },
				s => {
					// if s is a match, add it to the "concurrent queue"
					if (r.IsMatch (s)) {
						matches.Enqueue(s);
					}
				}
			);

			// convert list back to array now that we're done iterating
			string[] matches_str = matches.ToArray ();

			// check array for nulls. this was an issue before I was using concurrentqueue (back in the list days)
			// not needed anymore, but the code is left in JUST IN CASE :)
			// matches_str = matches_str.Where(x => !string.IsNullOrEmpty(x)).ToArray();

			// count number of matches and write friendly message
			int n_matches = matches_str.Length;
			if (noisy) {
				watch.Stop ();
				string time = watch.Elapsed.ToString(@"m\:ss");
				string tp = "Found " + n_matches + " matches in " + time + " minutes.";
				Console.WriteLine (OutputIndenter (tp, indent_lvl + 1));
			}
			return matches_str;
		}

	// KmerPositionSearch
	// finds position of kmer in multiple sequences
		public int[] KmerPositionSearch (string[] inputseqs, string kmer){
			// calculate length only once
			int nseqs = inputseqs.Length;
			// pre-allocate array for kmer positions
			int[] positions = new int[nseqs];
			// regex object for kmer
			for (int i = 0; i < nseqs; i++ ) {
				Match ri = Regex.Match (inputseqs [i], kmer);
				if (ri.Success) {
					int index = ri.Index;
					positions [i] = index;
				} else {
					positions [i] = -1;
				}
			}
			return positions;
		}
	
	// HitFilter
	// filters found sequences to exclude crappy sequences
		public string[] HitFilter (string[] inputseqs, string kmer, string contig, int overlap_cutoff = 0, bool noisy = true, int indent_lvl=0, string input_desc = "matching seqs"){
			// write helpful and descriptive message to user :)
			if (noisy) {
				string tp = "Filtering " + input_desc + ".";
				Console.WriteLine (OutputIndenter (tp, indent_lvl));
			}

			// get kmer positions of all inputseqs
			int[] kmerpos = KmerPositionSearch (inputseqs, kmer);

			// start up a bool so we know which sequences are shit; counters for explanatory output (WHY seqs were filtered)
			bool[] goodseqs = Enumerable.Repeat(true, inputseqs.Length).ToArray();
			int n_kmernotfound = 0;
			int n_toolong = 0;
			int n_tooshort = 0;
			int n_rightsidemismatch = 0;
			int n_noseqadd = 0;

			// go through each seq and check some stuff. checking from easy to hard
			for (int i = 0; i < kmerpos.Length; i++) {

				// check if kmer was not found in sequence
				if(kmerpos[i] < 0){
					goodseqs[i] = false;
					n_kmernotfound++;
					continue; // like "next i" in more sensible languages
				}

				// check if kmer is at BEGINNING of sequence, and nothing of value would be added
				if (kmerpos [i] == 0) {
					goodseqs[i] = false;
					n_noseqadd++;
					continue; 
				}

				// check if right side of sequence is longer than the current contig
				// this would produce errors for the next part
				if ((inputseqs [i].Length - kmerpos [i]) >= contig.Length) {
					goodseqs[i] = false;
					n_toolong++;
					continue;
				}

				// check if right side of sequence (including kmer) is long enough. important for contigs.
				string seq_right_i = inputseqs[i].Substring(kmer.Length + kmerpos[i]);
				if (seq_right_i.Length + kmer.Length < overlap_cutoff) {
					goodseqs[i] = false;
					n_tooshort++;
					continue;
				}

				// check if right side of sequence (after kmer) aligns with homologous part of contig
				string con_right_i = contig.Substring(kmer.Length, seq_right_i.Length);
				if (seq_right_i != con_right_i) {
					goodseqs[i] = false;
					n_rightsidemismatch++;
				}


			}

			// filter inputseqs by goodseqs
			// in R this would be "outputseqs <- inputseqs[goodseqs]", why is c# so verbose???
			string[] outputseqs = inputseqs.Where((x, index) => goodseqs[index]).ToArray();

			// filtering results messages
			if (noisy) {
				if (n_kmernotfound > 0) {
					string tp1 = "kmer was NOT found in " + n_kmernotfound.ToString () + " sequences. This is likely a bug.";
					Console.WriteLine (OutputIndenter (tp1, indent_lvl + 1));
				}
				if (n_noseqadd > 0) {
					string tp2 = "No new sequence added by " + n_noseqadd.ToString () + " sequences (kmer at position 0).";
					Console.WriteLine (OutputIndenter (tp2, indent_lvl + 1));
				}
				if (n_toolong > 0) {
					string tp6 = "Right side of kmer was too long in " + n_toolong.ToString () + " sequences.";
					Console.WriteLine (OutputIndenter (tp6, indent_lvl + 1));
				}
				if (n_tooshort > 0) {
					string tp7 = "Right side of kmer was shorter than cutoff in " + n_tooshort.ToString () + " sequences.";
					Console.WriteLine (OutputIndenter (tp7, indent_lvl + 1));
				}
				if (n_rightsidemismatch > 0) {
					string tp3 = "Right side of kmer did NOT align in " + n_rightsidemismatch.ToString () + " sequences.";
					Console.WriteLine (OutputIndenter (tp3, indent_lvl + 1));
				}
				int totalfiltered = goodseqs.Where (c => !c).Count ();
				int totalkept = goodseqs.Where (c => c).Count ();
				string tp4 = "Filtered " + totalfiltered + " sequences.";
				Console.WriteLine (OutputIndenter (tp4, indent_lvl + 1));
				string tp5 = "Kept: " + totalkept + " sequences.";
				Console.WriteLine (OutputIndenter (tp5, indent_lvl + 1));
			}

			// return sequences that passed filtering
			return outputseqs;
		
		}

	// Aligner
	// aligns an array of sequences using a shared kmer
		public string[] Aligner(string[] inputseqs, string kmer){
			// get kmer positions of all the seqs
			int[] kmerpos = KmerPositionSearch (inputseqs, kmer);
			int maxkmerpos = kmerpos.Max();

			// for each sequence, pad its left side and calculate total length
			int nseqs = inputseqs.Length;
			int[] seqlengths = new int[nseqs];
			for (int i = 0; i < nseqs; i++) {
				int pre_nchar = inputseqs [i].Length;
				// make a string of dashes equal to the max - this thing's start
				string lpad = new String ('-', maxkmerpos - kmerpos[i]);
				inputseqs [i] = lpad + inputseqs [i];
				seqlengths [i] = inputseqs [i].Length;
			}

			// for each sequence, pad the right side
			int maxlength = seqlengths.Max();
			for (int i = 0; i < nseqs; i++) {
				// make a string of dashes equal to maxlength - seqlengths[i]
				string rpad = new String ('-', maxlength - seqlengths[i]);
				inputseqs [i] = inputseqs [i] + rpad;
			}


			return(inputseqs);

		}

	// CountXinY
		// a dumb function that should not have to be written; faster than linq
		public int CountXinY(char x, char[] y){
			int count = 0;
			for (int i = 0; i < y.Length; i++) {
				if(y[i] == x){count++;}
			}
			return count;
		}

	// WhichMinMax
	// another dumbass function that's faster than using linq
	// !!!!! RETURNS -1 FOR A TIE !!!!!!
		public int WhichMinMax(int[] input, string minmax, bool tieOK = false){
			// find min or max values
			int output = 0; //the positional index of extremeval
			int extremeval = input [0]; // start extremeval at the first one

			for (int i = 0; i < input.Length; i++) {
				if (minmax == "max" && input [i] > extremeval) {
					extremeval = input [i];
					output = i;
				} else if (minmax == "min" && input [i] < output) {
					extremeval = input [i];
					output = i;
				}
			}
			// test if min or max value occurs multiple times (tie).
			if (tieOK == false) {
				int extremecount = 0;
				for (int i = 0; i < input.Length; i++) {
					if (input [i] == input[output]) {
						extremecount++;
					}
				}
				if (extremecount > 1) {
					output = -1;
				}
			}
			return output;
		}

	// MostCommonNT
	// From a column in an alignment, calculates the most common nucleotide.
	// returns '-' for bimodal columns, or all-gaps.
		public char MostCommonNT(char[] nts, bool nongap = true){
			// remove gaps from nts
			bool passthrough = false;
			if (nongap) {
				// check if everything is a gap
				if(nts.All(a => a == '-')){
					// set passthrough to true, skip rest of code
					passthrough = true;
				} else {
					// remake nts without gaps
					nts = nts.Where ((char a) => a != '-').ToArray();
				}
			}

			// set up output, initialize to failure
			// mcnt = "most common nucleotide"
			// initialize to gap so if passthrough is true (allgaps), you're done
			char mcnt = '-';

			if (passthrough == false) {
				// unique characters in array
				char[] elements = nts.Distinct ().ToArray ();
				// check if there's only one character in the array
				if (elements.Length == 1) {
					// if there's only one element in the whole column, that's the answer!
					mcnt = elements [0];
				} else {
					// find out counts of each element
					int[] counts = new int[elements.Length];
					for (int i = 0; i < elements.Length; i++) {
						counts [i] = CountXinY (elements [i], nts);
					}

					// get most common bp using counts, but check for bimodality too
					int mcnt_position = WhichMinMax (counts, "max");
					if (mcnt_position == -1) {
						// if there's a tie, output a gap, since there is no consensus.
						mcnt = '-';
					} else {
						mcnt = elements [mcnt_position];
					}
				}
			}

			return(mcnt);
		}

	// RowSums
	// computes row sums of a 2 dimentional array of ints
		public int[] RowSums(int[,] inmat){
			int ncol = inmat.GetLength(1);
			int nrow = inmat.GetLength(0);
			int[] output = new int[nrow];
			// for each row... 
			for (int i = 0; i < nrow; i++) {
				int sum = 0;
				for (int j = 0; j < ncol; j++) {
					sum = sum + inmat [i, j];
				}
				output [i] = sum;
			}
			return output;
		}

	// FilterAlignment
	// does some more filtering once sequences are aligned
		public string[] FilterAlignment(string[] seqs_aln, int fudge = 4, bool noisy = true, int indent_lvl=0, string input_desc = "aligned sequences"){

			// initialize matrix to store agreement values
			int nseqs = seqs_aln.Length;
			int seqlength = seqs_aln [1].Length;

			// write friendly messaage
			if (noisy) {
				string tp = "Filtering " +nseqs.ToString() + " " + input_desc + ".";
				Console.WriteLine (OutputIndenter (tp, indent_lvl));
			}


			//Matrix<int> agreementmat = Matrix<int>.Build.Dense (nseqs, seqlength);
			int[,] agreementmat = new int[nseqs, seqlength];

			// go through columns and rows, see if each value is the mode of its column
			for (int j = 0; j < seqlength; j++) {
				// get column of nucleotides at position j
				char[] nt_col_j = new char[nseqs];
				for (int i = 0; i < nseqs; i++) {
					nt_col_j [i] = seqs_aln [i] [j];
				}

				// find most common non-gap nucleotide for column j
				char mcnt = MostCommonNT(nt_col_j, nongap:true);

				// check which values at j are equal to the most common nucleotide
				for (int i = 0; i < nseqs; i++) {
					if (nt_col_j [i] == mcnt) {
						agreementmat [i, j] = 1;
					} else {
						agreementmat [i, j] = 0;
					}
				}
			}

			// compute agreement scores for each row
			int[] agreescores = RowSums(agreementmat);

			// calculate cutoff threshold; see which seqs are OK
			double threshold = agreescores.Average() - fudge;
			bool[] goodseqs = new bool[nseqs];
			for (int j = 0; j < nseqs; j++) {
				if (agreescores [j] > threshold) {
					goodseqs [j] = true;
				} else {
					goodseqs [j] = false;
				}
			}

			// tell user how many sequences we threw out
			if (noisy) {
				int totalfiltered = goodseqs.Where (c => !c).Count ();
				int totalkept = goodseqs.Where (c => c).Count ();
				string tp4 = "Filtered " + totalfiltered + " sequences.";
				Console.WriteLine (OutputIndenter (tp4, indent_lvl + 1));
				string tp5 = "Kept: " + totalkept + " sequences.";
				Console.WriteLine (OutputIndenter (tp5, indent_lvl + 1));
			}

			// throw out sequences that have fewer modal base-pairs than the average + fudge
			// again, no idea why or how this works
			string[] outputseqs = seqs_aln.Where((x, index) => goodseqs[index]).ToArray();
			return outputseqs;
		}
			
	// Consensus
		public string Consensus(string[] seqs_aln, string kmer, int coverage_threshold, int fudge=4, int indent_lvl=0, bool noisy = true, string input_desc = "aligned sequences"){
			int nseqs = seqs_aln.Length;
			// startup message
			if (noisy) {
				string tp = "Consensus of " + nseqs.ToString() + " " + input_desc + ".";
				Console.WriteLine (OutputIndenter (tp, indent_lvl));
			}

			// find starting position AND check for alignment
			// get kmer positions of all the seqs
			int[] kmerpos = KmerPositionSearch (seqs_aln, kmer);
			int startpos = new int();
			if (kmerpos.All (a => a == kmerpos [1])) {
				if (noisy) {
					string tp2 = "Confirmed: sequences have same kmer position.";
					Console.WriteLine (OutputIndenter (tp2, indent_lvl + 1));
				}
				startpos = kmerpos [0] - 1;

			} else {
				throw new Exception("ERROR - seqs input to consensus function are NOT aligned. This is a big bug.");
			}

			// move left from startpos and get consensus
			string output = "i";
			for (int p = startpos; p >= 0; p--) {
				// p = 'p'osition in alignment
				// get all nucleotides from alignment at p
				char[] nt_col_p = new char[nseqs];
				for (int i = 0; i < nseqs; i++) {
					nt_col_p [i] = seqs_aln [i] [p];
				}

				// figure out most comon nucleotide at p
				char nt_p = MostCommonNT(nt_col_p);

				// check if it's over the coverage threshold
				int cov_nt_p = nt_col_p.Where(x => x == nt_p).Count();
				if (cov_nt_p < coverage_threshold) {
					break;
				}

				// if gap, yer done. if not gap, add it to output
				if (nt_p != '-') {
					output = nt_p + output;
				} else {
					break;
				}
			}
			if (noisy) {
				string tp3 = "Made consensus of " + output.Length.ToString() + " nt to the left of kmer.";
				Console.WriteLine (OutputIndenter (tp3, indent_lvl + 1));
			}

			// remove the initialized i from the string
			output = output.Substring(0, output.Length-1);

			return output;
		}
	
	// ReverseComplement
		public string ReverseComplement(string inputseq, bool reverse=true){
			//Console.WriteLine ("raw_inputseq: " + inputseq);
			char[] in_ary = new char[inputseq.Length];
			in_ary = inputseq.ToCharArray ();
			char[] out_ary = new char[in_ary.Length];
			// complement first
			for (int i = 0; i < in_ary.Length; i++) {
				if (in_ary [i] == 'A') {
					out_ary [i] = 'T';
				} else if (in_ary [i] == 'T') {
					out_ary [i] = 'A';
				} else if (in_ary [i] == 'C') {
					out_ary [i] = 'G';
				} else if (in_ary [i] == 'G') {
					out_ary [i] = 'C';
				} else {
					out_ary [i] = '-';
				}
			}
			if (reverse == true) {
				out_ary = out_ary.Reverse().ToArray();
			}
			string outputseq = new string (out_ary);
			return(outputseq);
		}

	// Joke
	// tells a joke and tortures the user
	// in the next version, this will tell a joke at random.
	// todo: make a [3,n] string array, where the first index is question, response, answer. second index is joke - randomize that.
		public void Joke(){
			Console.WriteLine ("Q: why did the contig cross the road?");
			string answer = "not given";
			bool answered = false;
			while(answered == false){
				answer = Console.ReadLine ();
				if(answer.Equals("why", StringComparison.CurrentCultureIgnoreCase) || answer.Equals("why?", StringComparison.CurrentCultureIgnoreCase)){
					answered = true;
				} else { 
					Console.WriteLine("Invalid response. Expected: \"Why?\"");
				}
			}
			Console.WriteLine ("To get to the other side!");

			// give user time to read joke
			Console.ReadKey();
		}

	// CombineFwdRevSeqs
		public string[] CombineFwdRevSeqs(string[] fwd_seqs, string[] rev_seqs){
			// re-rc the rev seqs
			for (int i = 0; i < rev_seqs.Length; i++) {
				rev_seqs [i] = ReverseComplement (rev_seqs [i]);
			}
			// combine arrays
			string[] output = fwd_seqs.Concat(rev_seqs).ToArray();
			return output;
		}

	// GetBestContig
		public string GetBestContig(string[] filt_contigs, string kmer, bool include_kmer=false){
			int[] kmerpos = KmerPositionSearch (filt_contigs, kmer);
			int best_index = WhichMinMax (kmerpos, "max", tieOK: true);
			string best_contig = filt_contigs[best_index];
			if (include_kmer == false) {
				best_contig = best_contig.Substring (0, kmerpos [best_index]);
			}
			return best_contig;
		}

	}

}


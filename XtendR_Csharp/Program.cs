
// Contig extender (genome assembler) written by John L. Darcy
// Orinally written in R, this is the C# version. 
// email me at darcyj@colorado.edu 
// Sorry to everyone reading this source code who prefers curly braces on their own lines. SORRY.

using System;
using System.IO;
using System.Linq;
using CommandLine;


namespace XtendR_Csharp {

	// main code of the program - this is where the magic happens
	class MainClass {
		public static void Main (string[] args) {

			// get functions
			Functions functions = new Functions();

			// declare variables to be assigned to options, initialize to failure state for required.
			string raw_reads_path = "ERROR";
			string starting_contig_path = "ERROR";
			string output_path = "ERROR";
			int kmerlength = -1;
			int targetLength = -1;
			string contigs_path = "not used";
			int contig_overlap = 0;
			int consensus_threshold = -1;
			int min_kmer_hits = -1;
			int recovery_trim = -1;
			int nthreads = -1;
			bool quiet = false;
			bool citation = false;
			bool joke = false;
			bool help = false;

			// parse options
			var options = new ExtendR_Options();
			if (CommandLine.Parser.Default.ParseArguments (args, options)) {

				// re-assign values based on input
				raw_reads_path = options.raw_reads_path_in;
				starting_contig_path = options.starting_contig_path_in;
				kmerlength = options.kmerlength_in;
				output_path = options.output_path_in;
				targetLength = options.targetLength_in;
				contigs_path = options.contigs_path_in;
				contig_overlap = options.contig_overlap_in;
				consensus_threshold = options.consensus_threshold_in;
				min_kmer_hits = options.min_kmer_hits_in;
				recovery_trim = options.recovery_trim_in;
				nthreads = options.nthreads_in;
				quiet = options.quiet_in;
				citation = options.citation_in;
				joke = options.joke_in;
				help = options.help_in;

				// display help dialog
				if (help) {
					Console.Error.Write(options.getusage());
				}

				// tell a joke?
				if(joke){
					functions.Joke();
				}

				// handle my STUPID decition to use "noisy" in functions and "quiet" as input
				bool noisy = !quiet;

				// write parameters out so user can see
				if (noisy) {
					Console.WriteLine ("Parameter values and inputs:");
					Console.WriteLine (functions.OutputIndenter ("Raw reads: " + raw_reads_path, indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("Starting contig: " + starting_contig_path, indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("kmer length: " + kmerlength.ToString (), indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("Output path: " + output_path, indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("Target length: " + targetLength.ToString (), indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("Contigs path: " + contigs_path, indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("Consensus_threshold: " + consensus_threshold.ToString (), indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("Min kmer hits: " + min_kmer_hits.ToString (), indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("Recovery trim NT: " + recovery_trim.ToString (), indent_lvl: 1));
					Console.WriteLine (functions.OutputIndenter ("n threads: " + nthreads.ToString (), indent_lvl: 1));
				}

				// make sure output directory is set
				if (output_path == "ERROR") {
					Console.Error.Write("Output path not specified.");
				}

				// read in fasta sequences
				string[] seqs = functions.FastaReader(raw_reads_path, noisy, indent_lvl:0);

				// read in contigs, if there are some...?
				string[] contigs = new string[0];
				bool contigmode = false;

				if (contigs_path == "not used") {
					if (quiet == false) {
						Console.WriteLine (functions.OutputIndenter ("No contig filepath given - contig mode not used.", indent_lvl: 0));
					}
				} else {
					contigs = functions.FastaReader(contigs_path, noisy, indent_lvl:0);
					contigmode = true;
				}

				// read in starting contig
				string assembled_sequence = "";
				string[] starting_contig = functions.FastaReader(starting_contig_path, noisy, indent_lvl:0);
				if (starting_contig.Length != 1) {
					Console.Error.Write (functions.OutputIndenter ("Exactly 1 starting contig expected. Got " + starting_contig.Length.ToString() + ".", indent_lvl: 0));
				} else {
					assembled_sequence = starting_contig [0];
				}

				// do iterative assembly. 
				if (noisy) {
					Console.WriteLine ("Starting iterative assembly.");
				}
				int iteration = 0; 
				int addedlength = 0; // keep track of gains!
				int iter_indent = 1;
				int assembledlength = assembled_sequence.Length;

				while (addedlength < targetLength) {

					iteration++;

					// print starting iteration message
					if (noisy) {
						Console.WriteLine (functions.OutputIndenter ("Iteration " + iteration.ToString() + ".", iter_indent));
					}

					// set passthrough to false (this can be turned on to skip to the end of the loop where some messages are displayed
					bool passthrough = false;

					// initialize to_add, which is the new sequence to be added
					string to_add = "";

					// get kmer for current iteration
					string kmer = assembled_sequence.Substring(0, kmerlength);
					string rckmer = functions.ReverseComplement (kmer);

					// if contig mode, find matching contigs
					if (contigmode && passthrough==false) {
						// get and filter matching contigs
						string[] fwd_contigs = functions.ParallelKmerSearch(contigs, kmer, nthreads, noisy, iter_indent+1, input_desc:"forward contigs");
						string[] rev_contigs = functions.ParallelKmerSearch(contigs, rckmer, nthreads, noisy, iter_indent+1, input_desc:"reverse-complement contigs");
						string[] all_contigs = functions.CombineFwdRevSeqs (fwd_seqs: fwd_contigs, rev_seqs: rev_contigs);
						if (all_contigs.Length > 0) {
							string[] contigs_filtered = functions.HitFilter (all_contigs, kmer, assembled_sequence, overlap_cutoff: contig_overlap, noisy:noisy, indent_lvl:iter_indent+2, input_desc:"contigs");
							if (contigs_filtered.Length > 0) {

								to_add = functions.GetBestContig (contigs_filtered, kmer, include_kmer: false);
								passthrough = true;
							} else {
								Console.WriteLine (functions.OutputIndenter ("No useable contigs found this iteration.", iter_indent + 2));
								passthrough = false; //this is not necessary but i still like it.
							}
						} else {
							Console.WriteLine (functions.OutputIndenter ("No useable contigs found this iteration.", iter_indent + 2));
							passthrough = false; //this is not necessary but i still like it.
						}
					}

					// do fwd kmer search on raw sequences
					string[] fwd_matches = new string [0]; 
					if (passthrough == false) {
						fwd_matches = functions.ParallelKmerSearch (seqs, kmer, nthreads, noisy, iter_indent + 1);
						if (fwd_matches.Length > 0) {
							fwd_matches = functions.HitFilter(inputseqs:fwd_matches, kmer:kmer, contig:assembled_sequence, overlap_cutoff:0, noisy:noisy, indent_lvl:iter_indent+2);
						}
					}

					// did we get enough? if not, do rev kmer search
					string[] rev_matches = new string [0]; 
					if (passthrough == false && fwd_matches.Length < min_kmer_hits) {
						rev_matches = functions.ParallelKmerSearch (seqs, rckmer, nthreads, noisy, iter_indent + 1);
						if (rev_matches.Length > 0) {
							// reverse complement each so it has the fwd orientation

							// print all seqs for debug purposes
							//for (int i = 0; i < rev_matches.Length; i++) {
							//	Console.WriteLine ("i = " + i.ToString () + "; seq = " + rev_matches [i]);
							//}


							for (int i = 0; i < rev_matches.Length; i++) {
								rev_matches [i] = functions.ReverseComplement (inputseq:rev_matches [i], reverse: true);
							}
							rev_matches = functions.HitFilter(inputseqs:rev_matches, kmer:kmer, contig:assembled_sequence, overlap_cutoff:0, noisy:noisy, indent_lvl:iter_indent+2);
						}
					}

					// combine fwd and rev matches, get consensus
					if (passthrough == false) {
						string[] all_matches = fwd_matches.Concat(rev_matches).ToArray();
						//string[] all_matches = functions.CombineFwdRevSeqs (fwd_matches, rev_matches);
						if (all_matches.Length > consensus_threshold) {
							all_matches = functions.Aligner (all_matches, kmer);
							all_matches = functions.FilterAlignment (all_matches, 4, noisy, iter_indent + 1);
							if (all_matches.Length > consensus_threshold) {
								to_add = functions.Consensus (all_matches, kmer, consensus_threshold, 4, iter_indent + 1, noisy);
								passthrough = true;
							}
						}
					}

					// see if some sequence was actually added
					if (to_add.Length > 0) {
						assembled_sequence = to_add + assembled_sequence;
						// write message about how much was added and total length
						addedlength += to_add.Length;
						Console.WriteLine (functions.OutputIndenter ("Added " + to_add.Length.ToString() + " NTs, total length is now " + assembled_sequence.Length.ToString() + ", total added is " + addedlength.ToString(), iter_indent + 1));

					} else { 
						// write an error message, let the user know that we're trimming and trying again
						Console.WriteLine (functions.OutputIndenter ("No NTs were added this iteration. Trimming " + recovery_trim.ToString() + " NTs and trying again.", iter_indent + 1));
						assembled_sequence = assembled_sequence.Substring (recovery_trim, (assembled_sequence.Length - recovery_trim));
					}

					assembledlength = assembled_sequence.Length;
					string seqid = ">xtendrsequence_TotalLength=" + assembledlength.ToString() + "_AddedLength=" + addedlength.ToString();
					string[] fasta_out = {seqid, assembled_sequence};
					System.IO.File.WriteAllLines (path: output_path, contents: fasta_out);



				}

				Console.WriteLine ("Finished iterative assembly.");

				if (citation) {
					Console.WriteLine ("Citation: Darcy, JL. XtendR: Extending contigs in C#. 2017. Available at (github url goes here)");
				}
			
			}

		}
	}
}

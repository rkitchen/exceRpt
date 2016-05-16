package main;

import java.text.SimpleDateFormat;
import java.util.Calendar;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import transcriptome.CIGAR_2_PWM;
import fastqTools.FindAdapter;
import fastaTools.Fasta2Fastq;
import fastaTools.FastaHeaderGrep;
import fastqTools.Fastq2Fasta;
import fastqTools.FilterBySequenceLength;
import fastqTools.FilterFastxByHeaderList;
import fastqTools.GetSequenceLengths;
import fastqTools.MatchPairedEndSequences;
import fastqTools.ProcessFastqWithRandomBarcode;
import fastqTools.RemoveHomopolymers;
import footprintAlignments.ReadCoverage;

public class ExceRpt_Tools {

	public static final String VERSION = "1.2.3";
	
	public static final String OPT_PATH_DB_ANNOTATION = "A";
	public static final String OPT_PATH_DB_SAMPLE = "S";
	public static final String OPT_DB_TABLE_NAME = "t";
	public static final String OPT_PATH_INPUT = "i";
	public static final String OPT_PATH_OUTPUT = "o";
	public static final String OPT_PATH_SPECTRA_TANDEM = "i";
	public static final String OPT_PATH_SPECTRA_MZXML = "i";
	public static final String OPT_PATH_ANNOTATION = "a";
	public static final String OPT_FORMAT_GENCODE = "g";
	public static final String OPT_FORMAT_CUFFMERGE = "c";
	public static final String OPT_FORMAT_SWISSPROT = "s";
	public static final String OPT_FORMAT_TRINITY = "n";
	public static final String OPT_CHOICE_DB_FORCE_REFRESH = "F";

	/**
	 * Parse the command line arguments
	 * @param args
	 * @return
	 * @throws ParseException
	 */
	public static CommandLine parseArgs(String[] args, Options options) throws ParseException{
		CommandLineParser parser = new BasicParser();
		return parser.parse(options, args);	
	}



	public static void main(String[] args) throws Exception {
		String main = "NULL";
		if(args.length > 0){
			main = args[0].trim().toLowerCase(); 
		}

		if(main.equals("getsequencelengths")){
			GetSequenceLengths.main(args);
		}else if(main.equals("filtersequencesbylength")){
			FilterBySequenceLength.main(args);
		}else if(main.equals("removehomopolymerrepeats")){
			RemoveHomopolymers.main(args);
		}else if(main.equals("matchpairedendsequences")){
			MatchPairedEndSequences.main(args);
		}else if(main.equals("fastaheadergrep")){
			FastaHeaderGrep.main(args);
		}else if(main.equals("cigar_2_pwm")){
			CIGAR_2_PWM.main(args);
		}else if(main.equals("readcoverage")){
			ReadCoverage.main(args);
		}else if(main.equals("processendogenousalignments")){
			ProcessEndogenousAlignments.main(args);
		}else if(main.equals("quantifyendogenousalignments")){
			QuantifyEndogenousAlignments.main(args);
		}else if(main.equals("processexogenousalignments")){
			ProcessExogenousAlignments.main(args);
		}else if(main.equals("filterfastxbyidlist")){
			FilterFastxByHeaderList.main(args);
		}else if(main.equals("processfastqwithrandombarcode")){
			ProcessFastqWithRandomBarcode.main(args);
		}else if(main.equals("findadapter")){
			FindAdapter.main(args);
		}else if(main.equals("fastq2fasta")){
			Fastq2Fasta.main(args);
		}else if(main.equals("fasta2fastq")){
			Fasta2Fastq.main(args);
		}
		
		
		else{
			System.out.println("exceRpt_Tools version "+VERSION);
			System.out.println("");
			
			System.out.println("Usage:\t "+ExceRpt_Tools.ExceRpt_Tools_EXE_COMMAND+" <Command>");
			System.out.println("");

			System.out.println("Command: GetSequenceLengths            | Get the distribution of sequence lengths in a FASTA/Q file");
			System.out.println("         Fasta2Fastq                   | Convert FASTA sequence(s) to FASTQ sequence(s)");
			System.out.println("         Fastq2Fasta                   | Convert FASTQ sequence(s) to FASTA sequence(s)");
			System.out.println("         FastaHeaderGrep               | Filter fasta sequences based on the sequence ID");
			System.out.println("         FilterFastxByIDList           | Filter fasta/q sequences based on a list of sequence IDs");
			System.out.println("         FilterSequencesByLength       | Filter fasta or fastq sequences based on some maximum sequence length");
			System.out.println("         ProcessFastqWithRandomBarcode | Filter fasta or fastq sequences based on some maximum sequence length");
			System.out.println("         FindAdapter                   | Determine most likely 3' adapter sequence from fastq reads");
			System.out.println("         RemoveHomopolymerRepeats      | Filter fasta or fastq sequences based on sequence composition");
			System.out.println("         MatchPairedEndSequences       | Match paired-end fastq sequences based on readID");
			System.out.println("         CIGAR_2_PWM                   | Reads SAM/BAM alignments and converts the CIGAR strings to a position-weight matrix");
			System.out.println("         ReadCoverage                  | Reads SAM/BAM alignments to the TRANSCRIPTOME and calculates read coverage consistency");
			System.out.println("         ProcessEndogenousAlignments   | Process endogenous smallRNA alignments for the exceRpt pipeline");
			System.out.println("         QuantifyEndogenousAlignments  | Quantify endogenous smallRNA alignments for the exceRpt pipeline");
			System.out.println("         ProcessExogenousAlignments    | Process exogenous genomic alignments for the exceRpt pipeline");
			System.out.println();
		}


	}

	public static final String ExceRpt_Tools_EXE_COMMAND = "java -Xmx10G -jar exceRpt_Tools.jar";

	
	public static String getTime(){
		return((new SimpleDateFormat("yyyy/MM/dd HH:mm:ss")).format(Calendar.getInstance().getTime()));
	}
	
	public static void printOut(String message){ System.out.print(getTime()+" "+message); }
	public static void printErr(String message){ System.err.print(getTime()+" "+message); }
	public static void printLineOut(String message){ System.out.println(getTime()+" "+message); }
	public static void printLineErr(String message){ System.err.println(getTime()+" "+message); }
	
	public static void printProgressBar(int percent){
	    StringBuilder bar = new StringBuilder("[");
	    for(int i = 0; i < 50; i++){
	        if( i < (percent/2)){
	            bar.append("=");
	        }else if( i == (percent/2)){
	            bar.append(">");
	        }else{
	            bar.append(" ");
	        }
	    }
	    bar.append("]   " + percent + "%     ");
	    printErr(bar.toString()+"\r");
	}
	
}

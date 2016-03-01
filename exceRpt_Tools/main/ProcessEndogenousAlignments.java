package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeSet;

import net.sf.samtools.Cigar;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import sequenceTools.Sequence;


public class ProcessEndogenousAlignments {

	//private BufferedWriter _optimalAlignmentWriter;
	HashMap <String, Hairpin> miRNAprecursors = new HashMap <String, Hairpin>();


	/*public ProcessEndogenousAlignments(String pathForOptimalAlignments) throws IOException{
		_optimalAlignmentWriter = new BufferedWriter(new FileWriter(pathForOptimalAlignments+"/endogenousAlignments_Accepted.txt"));
	}*/
	public ProcessEndogenousAlignments(){}


	private ArrayList<String> _libraryPriority;
	private final String _libraryPriorityOptions = "miRNA, tRNA, piRNA, gencode, circRNA";
	private String _libraryPriorityString_default = "miRNA,tRNA,piRNA,gencode,circRNA";
	private String _libraryPriorityString = "";


	private boolean _readingGenomeMappedReads = true;


	/**
	 * Set the default library priority (miRNA,tRNA,piRNA,gencode,circRNA)
	 */
	private void setLibraryPriority(){
		setLibraryPriority(_libraryPriorityString_default);
	}



	/**
	 * Specify a new library priority
	 * @param priorityString new library string of the form miRNA,tRNA,piRNA,gencode,circRNA
	 */
	private void setLibraryPriority(String priorityString){
		_libraryPriority = new ArrayList<String>();
		_libraryPriorityString = "";
		String[] bits = priorityString.split(",");
		for(int i=0;i<bits.length;i++){
			String lib = validateLibrary(bits[i].trim());
			if(lib != null){
				_libraryPriority.add(lib);
				if(i > 0)
					_libraryPriorityString += ",";
				_libraryPriorityString += lib;
			}else
				ExceRpt_Tools.printLineErr("WARNING: Ignoring unrecognised library \'"+bits[i].trim()+"\'. Must be one of "+_libraryPriorityOptions);
		}

		if(_libraryPriorityString.equals(_libraryPriorityString_default)){
			ExceRpt_Tools.printLineErr("Using default library priorities: "+_libraryPriorityString_default.replace(","," > "));
		}else{
			ExceRpt_Tools.printLineErr("Overriding default library priorities: "+_libraryPriorityString.replace(","," > "));
		}

	}



	/**
	 *  checks input library string against valid options (case insensitive)
	 * @param lib library name to test
	 * @return properly formatted library name, or null if no match found
	 */
	private String validateLibrary(String lib){
		String validLibString = "";

		/* requires Java 1.7
		switch(lib.toLowerCase()){
		case "mirna": validLibString = "validLibString"; break;
		case "trna": validLibString = "tRNA"; break;
		case "pirna": validLibString = "piRNA"; break;
		case "gencode": validLibString = "gencode"; break;
		case "circrna": validLibString = "circRNA"; break;
		default: validLibString = null;
		}*/

		String libLower = lib.toLowerCase();
		if(libLower.equalsIgnoreCase("miRNA"))
			validLibString = "miRNA";
		else if(libLower.equalsIgnoreCase("tRNA"))
			validLibString = "tRNA";
		else if(libLower.equalsIgnoreCase("piRNA"))
			validLibString = "piRNA";
		else if(libLower.equalsIgnoreCase("gencode"))
			validLibString = "gencode";
		else if(libLower.equalsIgnoreCase("circRNA"))
			validLibString = "circRNA";

		return validLibString;
	}






	/*
	 * 
	 */
	/*public void tidyUp() throws IOException{
		_optimalAlignmentWriter.flush();
		_optimalAlignmentWriter.close();
	}*/




	/**
	 * 
	 * @param readsByLibrary
	 * @param refType
	 * @param alignment
	 * @return
	 */
	private HashMap<String, ArrayList<SAMRecord>> addReadAssignment(HashMap<String, ArrayList<SAMRecord>> readsByLibrary, String refType, SAMRecord alignment){
		if(!readsByLibrary.containsKey(refType))
			readsByLibrary.put(refType, new ArrayList<SAMRecord>());

		readsByLibrary.get(refType).add(alignment);
		return readsByLibrary;
	}




	private HashMap<String, Integer> _referenceID2index = new HashMap<String, Integer>();
	private int _indexCount = 1;

	/**
	 * Processes the current read:
	 *  1- parses the alignments 
	 *  2- prioritises the alignment based on miRNA>tRNA>piRNA>gencode>circularRNA
	 *  3- writes accepted alignment(s)
	 * @param thisRead
	 * @throws IOException 
	 */
	public void assignRead(HashMap<SAMRecord, String> thisRead) throws IOException{
		Iterator<Entry<SAMRecord,String>> it = thisRead.entrySet().iterator();
		Entry<SAMRecord,String> tmp;

		// store the alignment records and the libraries to which they map
		HashMap<String, ArrayList<SAMRecord>> readsByLibrary = new HashMap<String, ArrayList<SAMRecord>>();  
		ArrayList<SAMRecord> keepAlignments = new ArrayList<SAMRecord>();
		boolean readHasSenseAlignment = false;

		/*
		 * Loop through alignments and count libraries that capture alignments
		 */
		//int count = 0;
		while(it.hasNext()){
			tmp = it.next();
			boolean readNegativeStrand = tmp.getKey().getReadNegativeStrandFlag();
			if(!readNegativeStrand  &&  !readHasSenseAlignment)
				readHasSenseAlignment = true;

			// if this read is antisense, reverse complement it
			if(readNegativeStrand)
				tmp.getKey().setReadString(Sequence.reverseComplement(tmp.getKey().getReadString()));

			String thisLib = "balls";
			if(tmp.getValue().equals("miRNA")){
				if(!readNegativeStrand)
					thisLib = "miRNA_sense";
				else
					thisLib = "miRNA_antisense";
			}else if(tmp.getValue().equals("tRNA")){
				if(!readNegativeStrand)
					thisLib = "tRNA_sense";
				else
					thisLib = "tRNA_antisense";
			}else if(tmp.getValue().equals("piRNA")){
				if(!readNegativeStrand)
					thisLib = "piRNA_sense";
				else
					thisLib = "piRNA_antisense";
			}else if(tmp.getValue().equals("gencode")){
				if(!readNegativeStrand)
					thisLib = "gencode_sense";
				else
					thisLib = "gencode_antisense";
			}else if(tmp.getValue().equals("circRNA")){
				if(!readNegativeStrand)
					thisLib = "circRNA_sense";
				else
					thisLib = "circRNA_antisense";
			}

			
			readsByLibrary = addReadAssignment(readsByLibrary, thisLib, tmp.getKey());
			//System.out.println(tmp.getKey().getReadName()+"\ttmp.getValue()="+tmp.getValue()+"\t"+tmp.getKey().getReferenceName()+"\t"+tmp.getValue()+"\t"+(!readNegativeStrand)+"\t"+readHasSenseAlignment+"\t"+thisLib);
		}


		/*
		 * Logic block to prioritise alignments:
		 * - prefer sense alignments over antisense alignments
		 */
		String keptLibrary = "none";
		String thisLib = "";
		// Try sense alignments
		if(readHasSenseAlignment){
			Iterator<String> libraryIterator = _libraryPriority.iterator();
			//System.out.print("SENSE: ");
			while(libraryIterator.hasNext()){
				thisLib = libraryIterator.next();
				//System.out.print("\t"+thisLib);
				if(readsByLibrary.containsKey(thisLib+"_sense")){
					keepAlignments = readsByLibrary.get(thisLib+"_sense");
					keptLibrary = thisLib+"_sense";
					break;
				}
			}
		}else{  // Try antisense alignments
			Iterator<String> libraryIterator = _libraryPriority.iterator();
			//System.out.print("ANTISENSE: ");
			while(libraryIterator.hasNext()){
				thisLib = libraryIterator.next();
				//System.out.print("\t"+thisLib);
				if(readsByLibrary.containsKey(thisLib+"_antisense")){
					keepAlignments = readsByLibrary.get(thisLib+"_antisense");
					keptLibrary = thisLib+"_antisense";
					break;
				}
			}
		}
		//System.out.println("\tkeptLib: "+keptLibrary);


		/*
		 * If this is a miRNA aligned read, check to see if it overlaps any known mature sequences
		 */
		Iterator<SAMRecord> it2;
		SAMRecord tmp2;
		if(keptLibrary.startsWith("miRNA_")){
			ArrayList<SAMRecord> precursorAlignments = new ArrayList<SAMRecord>();
			//int index = 0;
			it2 = keepAlignments.iterator();

			// if this is a miRNA alignment, run through all alignments and check for a hit to a mature miRNA 
			while(it2.hasNext()){
				//index ++;
				tmp2 = it2.next();
				String[] refIDbits = tmp2.getReferenceName().split(":");
				
				String mapsTo = refIDbits[1];
				int startAt = 2;
				if(_forceLibrary != null){	
					mapsTo = refIDbits[0];
					startAt = 1;
				}
				for(int i=startAt;i<refIDbits.length;i++)
					mapsTo = mapsTo.concat(":"+refIDbits[i]);

				// if this read maps to a mature miRNA, 
				if(miRNAprecursors.get(mapsTo).getMatureOverlapForRead(tmp2.getReadName(), tmp2.getAlignmentStart(), tmp2.getAlignmentEnd(), keptLibrary.endsWith("_antisense")).length() > 0){
					if(keptLibrary.endsWith("_sense"))
						keptLibrary = "miRNAmature_sense";
					else
						keptLibrary = "miRNAmature_antisense";
					break;
				}else{
					precursorAlignments.add(tmp2);
				}
			}

			if(keptLibrary.startsWith("miRNAmature_")){ // if there are mature alignments, loop through them again and remove non-mature alignments
				Iterator<SAMRecord> it3 = precursorAlignments.iterator();
				while(it3.hasNext()){
					keepAlignments.remove(it3.next());
				}

			}else if(keptLibrary.startsWith("miRNA_")){ // if no mature sequence explains this read, explicitly assign it to the precursor library
				if(keptLibrary.endsWith("_sense"))
					keptLibrary = "miRNAprecursor_sense";
				else
					keptLibrary = "miRNAprecursor_antisense";
			}
		}


		/*
		 * 
		 */
		it2 = keepAlignments.iterator();
		//String lastLib = "";
		int counter = 0;

		//String referenceIDs = "";
		TreeSet<String> refIDs = new TreeSet<String>();

		while(it2.hasNext()){
			counter ++;
			tmp2 = it2.next();
			//System.out.println(tmp2.getReferenceName());
			//String insertSequence = tmp2.getReadString();


			/*
			 * Parse the SAM reference string:
			 *   [genome|nogenome]:[library]:[referenceID]
			 * e.g.
			 *   genome:gencode:ENST00000408108.1:miRNA:MIR486-201:GN=ENSG00000221035.1 
			 */
			String[] refIDbits = tmp2.getReferenceName().split(":");
			
			String mapsTo = refIDbits[1];
			int startAt = 2;
			if(_forceLibrary != null){	
				mapsTo = refIDbits[0];
				startAt = 1;
			}
			for(int i=startAt;i<refIDbits.length;i++)
				mapsTo = mapsTo.concat(":"+refIDbits[i]);

			String isGenomeMapped = "genome";
			if(!_readingGenomeMappedReads)
				isGenomeMapped = "nogenome";

			/*
			 *  For the FIRST alignment, print the read info:
			 */
			if(counter == 1)
				System.out.print(tmp2.getReadName()+"\t"+tmp2.getReadString()+"\t"+keptLibrary+"\t"+isGenomeMapped+"\t");
			//else
			//System.out.print(" | ");
			//referenceIDs = referenceIDs.concat(" | ");


			/*
			 * If this is a miRNA aligned read that HAS an overlap with a mature sequence, change the library and reference it to that of the mature
			 * 
			 * -- Suppress this stuff, save allocating reads for counting until we read the accepted alignments
			 */
			//boolean writtenThisAlignment = false;
			String matureID = "NA";
			if(keptLibrary.startsWith("miRNAmature_")){
				matureID = miRNAprecursors.get(mapsTo).getMatureOverlapForRead(tmp2.getReadName(), tmp2.getAlignmentStart(), tmp2.getAlignmentEnd(), keptLibrary.endsWith("_antisense"));
				if(matureID.length() > 0){
					mapsTo = matureID;
				}
			}

			// if this is a new referenceID, add it to the sorted list:
			if(!refIDs.contains(mapsTo)){
				refIDs.add(mapsTo);
			}

		}

		// create reference key
		String referenceIDs = "";
		boolean first = true;
		for(String id:refIDs){
			if(first)
				referenceIDs = referenceIDs.concat(id);
			else
				referenceIDs = referenceIDs.concat(" | "+id);
			first = false;
		}

		if(!_referenceID2index.containsKey(referenceIDs)){
			_referenceID2index.put(referenceIDs, _indexCount);
			_dictionaryWriter.write(_indexCount+"\t"+referenceIDs+"\n");

			System.out.println(_indexCount);
			_indexCount ++;
		}else{
			System.out.println(_referenceID2index.get(referenceIDs));
		}


	}




	private BufferedWriter _dictionaryWriter; 

	/**
	 * Reads the read alignments
	 * reference ID format is
	 *   [genome|nogenome]:[library]:[referenceID]
	 * e.g.
	 *   genome:miRNA:hsa-mir-486-1:MI0002470:Homo:sapiens:miR-486:stem-loop
	 * 
	 * @param path_readAlignments
	 * @throws IOException 
	 */

	public boolean read_Reads(SAMFileReader inputSam) throws IOException{
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		//inputSam.setValidationStringency(ValidationStringency.LENIENT);

		// TODO: put this if/else back when we're done testing
		//if(inputSam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.queryname)){
		SAMRecord thisRecord;
		HashMap<SAMRecord, String> thisRead = new HashMap<SAMRecord, String>(); 
		SAMRecordIterator it = inputSam.iterator();
		String lastReadID = null;
		while(it.hasNext()){
			//count++;
			thisRecord = it.next();

			//System.out.println(thisRecord.getReferenceName());
			if(!thisRecord.getReadName().equals(lastReadID)  &&  lastReadID != null){
				// new, non first
				assignRead(thisRead);
				thisRead = new HashMap<SAMRecord, String>();
			}

			// put the SAM record into the map with the library type as the value
			//thisRead.put(thisRecord, thisRecord.getReferenceName().split(":")[1]);
			
			if(_forceLibrary != null)
				thisRead.put(thisRecord, _forceLibrary);
			else
				thisRead.put(thisRecord, thisRecord.getReferenceName().split(":")[0]);
			lastReadID = thisRecord.getReadName();
			
		}
		// assign the final read!
		assignRead(thisRead);

		/*}else{
		Thunder.printLineErr("ERROR: Input SAM file must be sorted by readID");
		inputSam.close();
		return false;
		}*/

		inputSam.close();
		return true;
	}


	public boolean read_Reads(File path_readAlignments_genome, File path_readAlignments_transcriptome, String outputPath) throws IOException{

		_dictionaryWriter = new BufferedWriter(new FileWriter(new File(outputPath)));

		// Read the alignments for genome and transcriptome mapped reads
		if(path_readAlignments_genome != null){
			ExceRpt_Tools.printLineErr("Processing transcriptome alignments for genome-mapped reads...");
			_readingGenomeMappedReads = true;
			read_Reads(new SAMFileReader(path_readAlignments_genome));
			_dictionaryWriter.flush();
		}

		// Read the alignments for transcriptome mapped reads
		ExceRpt_Tools.printLineErr("Processing transcriptome alignments for genome-UNmapped reads...");
		_readingGenomeMappedReads = false;
		read_Reads(new SAMFileReader(path_readAlignments_transcriptome));
		_dictionaryWriter.flush();

		_dictionaryWriter.close();
		return true;
	}


	/**
	 * 
	 * @param path_hairpin2genome
	 * @param path_mature2hairpin
	 */
	public void read_miRNAinfo(File path_hairpin2genome, File path_mature2hairpin){
		/*
		 * Read the hairpin alignments to the genome
		 */
		SAMFileReader inputSam = new SAMFileReader(path_hairpin2genome);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMRecord thisRecord;
		SAMRecordIterator it = inputSam.iterator();
		while(it.hasNext()){
			thisRecord = it.next();
			//System.out.println(thisRecord.getReferenceName()+"\t"+thisRecord.getReadName()+"\t"+thisRecord.getReadString());
			if(!miRNAprecursors.containsKey(thisRecord.getReadName()))
				miRNAprecursors.put(thisRecord.getReadName(), new Hairpin(thisRecord.getReadName()));

			if(!thisRecord.getReadUnmappedFlag())
				miRNAprecursors.get(thisRecord.getReadName()).addGenomicAlignment(new Alignment(thisRecord.getReferenceName(), thisRecord.getAlignmentStart(), thisRecord.getAlignmentEnd(), thisRecord.getReadNegativeStrandFlag(), thisRecord.getCigar(), thisRecord.getReadString()));
		}
		inputSam.close();

		/*Iterator<String> it2 = miRNAprecursors.keySet().iterator();
		while(it2.hasNext()){
			Hairpin tmp = miRNAprecursors.get(it2.next());
			System.out.println(tmp.getID()+"\t"+tmp.getNumberOfGenomicAlignments());
		}*/


		/*
		 * Read the mature sequences
		 */
		inputSam = new SAMFileReader(path_mature2hairpin);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		it = inputSam.iterator();
		while(it.hasNext()){
			thisRecord = it.next();

			//System.out.println(thisRecord.getReferenceName()+"\t"+thisRecord.getReadName());

			if(!thisRecord.getReadUnmappedFlag()  &&  miRNAprecursors.containsKey(thisRecord.getReferenceName()))
				miRNAprecursors.get(thisRecord.getReferenceName()).addMatureMiRNA(new Alignment(thisRecord.getReadName(), thisRecord.getAlignmentStart(), thisRecord.getAlignmentEnd(), thisRecord.getReadNegativeStrandFlag(), thisRecord.getCigar(), thisRecord.getReadString()));


		}
		inputSam.close();
	}

	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName(".SAM or .BAM").hasArg().withDescription("Path to HAIRPIN alignments to the genome").create("hairpin2genome"));
		options.addOption(OptionBuilder.withArgName(".SAM or .BAM").hasArg().withDescription("Path to MATURE alignments to the hairpins").create("mature2hairpin"));
		options.addOption(OptionBuilder.withArgName(".SAM or .BAM").hasArg().withDescription("Path to read alignments to the genome AND transcriptome [MUST be sorted by ReadID]").create("genomeMappedReads"));
		options.addOption(OptionBuilder.withArgName(".SAM or .BAM").hasArg().withDescription("Path to read alignments to the transcriptome ONLY [MUST be sorted by ReadID]").create("transcriptomeMappedReads"));
		//options.addOption(OptionBuilder.withArgName(".stats").hasArg().withDescription("[optional] Path to random barcode stats").create("randombarcode"));
		options.addOption(OptionBuilder.withArgName("csv list").hasArg().withDescription("[optional] Library priorities for quantification. Comma separated list of libraries in *descending* order of importance. Default: miRNA,tRNA,piRNA,gencode,circRNA - this can also be used to suppress libraries during quantification.").create("libPriority"));
		options.addOption(OptionBuilder.withArgName("path").hasArg().withDescription("File to write the alignment dictionary into").create("dict"));
		options.addOption(OptionBuilder.withArgName("string").hasArg().withDescription("[optional] Force library for all alignments.").create("forceLib"));
		return options;
	}


	public static void main(String[] args) throws ParseException, IOException {
		/*String hairpin2genome = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MICRO_RNA/miRBase_v21_hairpin_hsa_hg19_aligned.sam";
		String mature2hairpin = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MICRO_RNA/miRBase_v21_mature_hairpin_hsa_aligned.sam";
		String readsPath_GnT = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_genomeMapped_transcriptome_Aligned.out.bam";
		String readsPath_T = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_genomeUnmapped_transcriptome_Aligned.out.bam";
		String output_dictionary = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/alignments.dict";
		args = new String[]{"ProcessEndogenousAlignments",
				"--hairpin2genome",hairpin2genome,
				"--mature2hairpin",mature2hairpin,
				"--genomeMappedReads",readsPath_GnT,
				"--transcriptomeMappedReads",readsPath_T,
				"--dict",output_dictionary,
				"--libPriority","miRNA,tRNA,piRNA,gencode,circRNA"
		};*/

		/*String hairpin2genome = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MICRO_RNA/miRBase_v21_hairpin_hsa_hg19_aligned.sam";
		String mature2hairpin = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MICRO_RNA/miRBase_v21_mature_hairpin_hsa_aligned.sam";
		String readsPath_T = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/TEST.sam";
		String output_dictionary = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/alignments.dict";
		args = new String[]{"ProcessEndogenousAlignments",
				"--hairpin2genome",hairpin2genome,
				"--mature2hairpin",mature2hairpin,
				"--transcriptomeMappedReads",readsPath_T,
				"--dict",output_dictionary,
				"--libPriority","miRNA,tRNA,piRNA,gencode,circRNA"
		};*/
		
		/*String hairpin2genome = "/Users/robk/Downloads/miRNA_precursor2genome.sam";
		String mature2hairpin = "/Users/robk/Downloads/miRNA_mature2precursor.sam";
		String readsPath_T = "/Users/robk/Downloads/exogenous_miRBase_Aligned.out.bam";
		String output_dictionary = "/Users/robk/Downloads/alignments.dict";
		args = new String[]{"ProcessEndogenousAlignments",
				"--hairpin2genome",hairpin2genome,
				"--mature2hairpin",mature2hairpin,
				"--transcriptomeMappedReads",readsPath_T,
				"--dict",output_dictionary,
				"--forceLib","miRNA"
		};*/

		CommandLine cmdArgs = ExceRpt_Tools.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption("hairpin2genome") && cmdArgs.hasOption("mature2hairpin") && cmdArgs.hasOption("transcriptomeMappedReads") && cmdArgs.hasOption("dict")){
			//if(cmdArgs.hasOption("hairpin2genome") && cmdArgs.hasOption("mature2hairpin") && cmdArgs.hasOption("reads2all")){
			//ProcessEndogenousAlignments engine = new ProcessEndogenousAlignments(new File(cmdArgs.getOptionValue("hairpin2genome")), new File(cmdArgs.getOptionValue("mature2hairpin")));
			//ProcessEndogenousAlignments engine = new ProcessEndogenousAlignments(cmdArgs.getOptionValue("outputPath"));
			ProcessEndogenousAlignments engine = new ProcessEndogenousAlignments();

			// Read barcode stats if they exist
			//if(cmdArgs.hasOption("randombarcode")){
			//	Thunder.printLineErr("Reading random barcode stats");
			//	engine.readRandomBarcodeInfo(cmdArgs.getOptionValue("randombarcode"));
			//}

			// Check for library priority change
			if(cmdArgs.hasOption("libPriority")){
				engine.setLibraryPriority(cmdArgs.getOptionValue("libPriority"));
				//Thunder.printLineErr("Overriding default library priorities: "+engine._libraryPriorityString);
			}else{
				engine.setLibraryPriority();
				//Thunder.printLineErr("Using default library priorities: "+engine._libraryPriorityString);
			}

			if(cmdArgs.hasOption("forceLib"))
				engine.setForceLibrary(cmdArgs.getOptionValue("forceLib"));
			
			// Read hairpin alignments to the genome and mature alignments to the hairpins
			ExceRpt_Tools.printLineErr("Reading miRNA annotation info");
			engine.read_miRNAinfo(new File(cmdArgs.getOptionValue("hairpin2genome")), new File(cmdArgs.getOptionValue("mature2hairpin")));

			// read and process the RNA-seq read alignments
			ExceRpt_Tools.printLineErr("Processing RNA-seq alignments");
			File genomeMapped = null;
			if(cmdArgs.hasOption("genomeMappedReads"))
				genomeMapped = new File(cmdArgs.getOptionValue("genomeMappedReads"));
			boolean cont = engine.read_Reads(genomeMapped, new File(cmdArgs.getOptionValue("transcriptomeMappedReads")), cmdArgs.getOptionValue("dict"));

			if(cont)
				ExceRpt_Tools.printLineErr("Done!");
			else
				ExceRpt_Tools.printLineErr("Something went wrong...");
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(ExceRpt_Tools.ExceRpt_Tools_EXE_COMMAND+" ProcessEndogenousAlignments", getCmdLineOptions());
			System.err.println();
		}
	}
	
	private String _forceLibrary = null;
	public void setForceLibrary(String lib){ _forceLibrary = lib; }
}



class Hairpin{
	private String _id; 
	public Hairpin(String id){
		_id = id;
	}

	public ArrayList<Alignment> _genomicAlignments = new ArrayList<Alignment>();
	public HashMap<String, Alignment> _matureMiRNAs = new HashMap<String, Alignment>();
	//public HashMap<String, Alignment> _reads = new HashMap<String, Alignment>();

	//public void addGenomicAlignment(String chr, int start, int stop, boolean isNegativeStrand){
	public void addGenomicAlignment(Alignment newAlignment){ _genomicAlignments.add(newAlignment); }
	public void addMatureMiRNA(Alignment newAlignment){
		if(_matureMiRNAs.containsKey(newAlignment.getReferenceID())){
			if(_matureMiRNAs.get(newAlignment.getReferenceID()).isNegativeStrand()  &&  !newAlignment.isNegativeStrand())
				_matureMiRNAs.put(newAlignment.getReferenceID(), newAlignment); //if the existing alignment is antisense to the hairpin, replace it with a sense alignment				
		}else{
			_matureMiRNAs.put(newAlignment.getReferenceID(), newAlignment);
		}

	}
	//public void addRead(String readID, int start, int stop, boolean isNegativeStrand){ 
	//	_reads.put(readID, new Alignment(_id, start, stop, isNegativeStrand)); 
	//}


	public String getID(){ return _id; }
	public int getNumberOfGenomicAlignments(){ return _genomicAlignments.size(); }
	public ArrayList<Alignment> getGenomicAlignments(){ return _genomicAlignments; }


	public String getMatureOverlapForRead(String readID, int start, int stop, boolean isNegativeStrand){
		double minFractionForOverlap = 0.8;
		String result = "";
		String sep = "";
		Iterator<String> matures = _matureMiRNAs.keySet().iterator();
		//System.out.println("\n>>"+readID+"\t"+_reads.get(readID).getReferenceID()+"\t"+_reads.get(readID).getStart()+"\t"+_reads.get(readID).getStop()+"\t"+_reads.get(readID).isNegativeStrand());
		while(matures.hasNext()){
			String tmpID = matures.next();
			_matureMiRNAs.get(tmpID).setReferenceID(_id);
			//System.out.println(">"+tmpID+"\t"+_matureMiRNAs.get(tmpID).getReferenceID()+"\t"+_matureMiRNAs.get(tmpID).getStart()+"\t"+_matureMiRNAs.get(tmpID).getStop()+"\t"+_matureMiRNAs.get(tmpID).isNegativeStrand());

			if((new Alignment(_id, start, stop, isNegativeStrand)).overlaps(_matureMiRNAs.get(tmpID), minFractionForOverlap)){
				result = result.concat(sep+tmpID);
				sep = "|";
			}
		}
		return result;
	}

}

class Alignment{
	private String _referenceID, _querySequence;
	private int _start, _stop;
	private boolean _isNegativeStrand;
	private Cigar _cigar;
	public Alignment(String referenceID, int start, int stop, boolean isNegativeStrand, Cigar cigar, String querySequence){
		_referenceID = referenceID;
		_start = start;
		_stop = stop;
		_isNegativeStrand = isNegativeStrand;
		_cigar = cigar;
		_querySequence = querySequence;
	}
	public Alignment(String referenceID, int start, int stop, boolean isNegativeStrand){
		_referenceID = referenceID;
		_start = start;
		_stop = stop;
		_isNegativeStrand = isNegativeStrand;
	}

	public void setReferenceID(String newID){ _referenceID = newID; }
	public String getReferenceID(){ return _referenceID; }
	public int getStart(){ return _start; }
	public int getStop(){ return _stop; }
	public boolean isNegativeStrand(){ return _isNegativeStrand; }
	public Cigar getCIGAR(){ return _cigar; }
	public String getSequence(){ return _querySequence; }

	public boolean overlaps(Alignment otherAlignment, double minFracOverlap){
		boolean result = false;
		if(_referenceID.equals(otherAlignment.getReferenceID())){
			if(_isNegativeStrand == otherAlignment.isNegativeStrand()){
				if((_start >= otherAlignment.getStart()  &&  _start <= otherAlignment.getStop())  ||  (_stop >= otherAlignment.getStart()  &&  _stop <= otherAlignment.getStop())){
					int diff = Math.abs(_start - otherAlignment.getStart()) + Math.abs(_stop - otherAlignment.getStop()); 
					int max = Math.max(_stop, otherAlignment.getStop()) - Math.min(_start, otherAlignment.getStart());
					//System.out.print(1.0 - (diff / (max+0.0))+"\t"+minFracOverlap+"\t"+(1.0 - (diff / (max+0.0)) >= minFracOverlap));
					if((1.0 - (diff / (max+0.0)) >= minFracOverlap))
						result = true;
				}
			}
		}
		return result;
	}
}




class ValueComparator_DoubleMap implements Comparator<String> {

	Map<String, Double> base;
	public ValueComparator_DoubleMap(Map<String, Double> base) {
		this.base = base;
	}

	// Note: this comparator imposes orderings that are inconsistent with equals.    
	public int compare(String a, String b) {
		if (base.get(a) >= base.get(b)) {
			return -1;
		} else {
			return 1;
		} // returning 0 would merge keys
	}
}


/**
 * Compares values based on the third entry of an expression count array:
 * {uniqueReadCount, totalReadCount, multimapAdjustedReadCount, multimapAdjustedBarcodeCount}
 * @author robk
 *
 */
class Comparator_SortByAdjustedReadCount implements Comparator<String> {

	Map<String, double[]> base;
	public Comparator_SortByAdjustedReadCount(Map<String, double[]> base) {
		this.base = base;
	}

	// Note: this comparator imposes orderings that are inconsistent with equals.    
	public int compare(String a, String b) {
		if (base.get(a)[2] >= base.get(b)[2]) {
			return -1;
		} else {
			return 1;
		} // returning 0 would merge keys
	}
}


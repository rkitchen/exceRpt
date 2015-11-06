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
import java.util.TreeMap;

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

/**
 * Holds random barcode info for each unique insert sequence.
 * Also used to hold the accepted alignment info for these reads   
 * @author robk
 *
 */
class Insert{
	private String _insertSequence;
	private int _nInsertSequences;
	//private int _nUniqueBarcodes_5p, _nUniqueBarcodes_3p, _nUniqueBarcodes_5p3p;
	//private double _KLdivergence_5p, _KLdivergence_3p, _KLdivergence_5p3p;

	Insert(String insertSequence, int nInsertSequences){
		_insertSequence = insertSequence;
		_nInsertSequences = nInsertSequences;
	}
	/*void addRandomBarcode5p(int nUniqueBarcodes_5p, double KLdivergence_5p){
		_nUniqueBarcodes_5p = nUniqueBarcodes_5p;
		_KLdivergence_5p = KLdivergence_5p;
		_has5p = true;
	}
	void addRandomBarcode3p(int nUniqueBarcodes_3p, double KLdivergence_3p){
		_nUniqueBarcodes_3p = nUniqueBarcodes_3p;
		_KLdivergence_3p = KLdivergence_3p;
		_has3p = true;
	}
	void addRandomBarcode5p3p(int nUniqueBarcodes_5p3p, double KLdivergence_5p3p){
		_nUniqueBarcodes_5p3p = nUniqueBarcodes_5p3p;
		_KLdivergence_5p3p = KLdivergence_5p3p;
		_has5p3p = true;
	}*/

	private HashMap <String, Integer> _barcodeCounts = new HashMap <String, Integer>(); 
	void addRead(String readID){ 
		_nInsertSequences++;

		String barcode_5p="", barcode_3p="";//, barcode_5p3p="";
		String barcode = "";
		// Check for 5' barcode
		String[] bits = readID.split("\\{5p");
		if(bits.length == 2){
			barcode_5p = bits[1].split("\\}")[0];
			barcode = barcode_5p;
			_has5p = true;
		}
		// Check for 3' barcode
		bits = readID.split("\\{3p");
		if(bits.length == 2){
			barcode_3p = bits[1].split("\\}")[0];
			barcode += barcode_3p;
			_has3p = true;
		}
		// Check for 5' AND 3' barcodes
		if(_has5p && _has3p){
			_has5p3p = true;
			//barcode_5p3p = barcode_5p+barcode_3p;
		}

		// add this random barcode to the list
		if(_has5p || _has3p){
			if(!_barcodeCounts.containsKey(barcode))
				_barcodeCounts.put(barcode, 1);
			else
				_barcodeCounts.put(barcode, _barcodeCounts.get(barcode)+1);
		}

	}

	public int getNumberOfRandomBarcodes(){
		return _barcodeCounts.size();
	}


	private boolean _has5p=false, _has3p=false, _has5p3p=false;
	boolean hasRandomBarcode5p(){ return _has5p; }
	boolean hasRandomBarcode3p(){ return _has3p; }
	boolean hasRandomBarcode5p3p(){ return _has5p3p; }

	String getInsertSequence(){ return _insertSequence; }
	int getInsertSequenceCount(){ return _nInsertSequences; }
	//int getBarcodeCount_5p(){ return _nUniqueBarcodes_5p; }
	//int getBarcodeCount_3p(){ return _nUniqueBarcodes_3p; }
	//int getBarcodeCount_5p3p(){ return _nUniqueBarcodes_5p3p; }
	//double getKL_5p(){ return _KLdivergence_5p; }
	//double getKL_3p(){ return _KLdivergence_3p; }
	//double getKL_5p3p(){ return _KLdivergence_5p3p; }

	int getBarcodeCount(){
		/*if(_has5p3p)
			return _nUniqueBarcodes_5p3p;
		else if(_has5p)
			return _nUniqueBarcodes_5p;
		else if(_has3p)
			return _nUniqueBarcodes_3p;
		else
			return 0;*/
		return _barcodeCounts.size();
	}
	/*double getKL(){
		if(_has5p3p)
			return _KLdivergence_5p3p;
		else if(_has5p)
			return _KLdivergence_5p;
		else if(_has3p)
			return _KLdivergence_3p;
		else
			return 0;
	}*/

	//
	//
	//

	private String _mapsToLib;
	void setMapsToLib(String library){ _mapsToLib = library; }
	String mapsToLib(){ return _mapsToLib; }

	private ArrayList<String> _mapsToReferences = new ArrayList<String>();
	void addReferenceAlignment(String referenceID){
		if(!_mapsToReferences.contains(referenceID))
			_mapsToReferences.add(referenceID);
	}
	ArrayList<String> getReferenceAlignments(){ return _mapsToReferences; }

}

public class ProcessEndogenousAlignments {

	private BufferedWriter _optimalAlignmentWriter;
	HashMap <String, Hairpin> miRNAprecursors = new HashMap <String, Hairpin>();


	public ProcessEndogenousAlignments(String pathForOptimalAlignments) throws IOException{
		_optimalAlignmentWriter = new BufferedWriter(new FileWriter(pathForOptimalAlignments+"/endogenousAlignments_Accepted.txt"));
	}

	//public ProcessEndogenousAlignments(String pathForOptimalAlignments, String pathToRandomBarcodeStats) throws IOException{
	//	_optimalAlignmentWriter = new BufferedWriter(new FileWriter(pathForOptimalAlignments+"/endogenousAlignments_Accepted.txt"));
	//	readRandomBarcodeInfo(pathToRandomBarcodeStats);
	//}

	private boolean _hasBarcodeStats = false;
	private HashMap<String, Insert> _allInserts = new HashMap<String, Insert>(); 

	
	

	private ArrayList<String> _libraryPriority;
	private final String _libraryPriorityOptions = "miRNA, tRNA, piRNA, gencode, circRNA";
	private String _libraryPriorityString_default = "miRNA,tRNA,piRNA,gencode,circRNA";
	private String _libraryPriorityString = "";
	
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
	


	/**
	 * Read summary info from the random barcodes (if they are present)
	 * 
	 * @param path
	 * @throws IOException
	 */
	/*private void readRandomBarcodeInfo(String path) throws IOException{
		_hasBarcodeStats = true;

		BufferedReader barcodeStatsReader = new BufferedReader(new FileReader(path));
		String line;
		String[] bits, header;

		// Look at stats file header row to determine the arrangement of barcodes (5' and/or 3')
		header = barcodeStatsReader.readLine().split("\t");
		int index_5p=-1, index_3p=-1, index_5p3p=-1;
		for(int i=0;i<header.length;i++){
			if(header[i].equalsIgnoreCase("nUniqueBarcodes_5p"))
				index_5p = i;
			else if(header[i].equalsIgnoreCase("nUniqueBarcodes_3p"))
				index_3p = i;
			else if(header[i].equalsIgnoreCase("nUniqueBarcodes_5p3p"))
				index_5p3p = i;
		}

		// read the rest of the stats file
		while((line=barcodeStatsReader.readLine()) != null){
			// InsertSequence	nInsertSequences	nUniqueBarcodes_5p	KLdivergence_5p	topBarcode_5p	secondBarcode_5p	thirdBarcode_5p	nUniqueBarcodes_3p	KLdivergence_3p	topBarcode_3p	secondBarcode_3p	thirdBarcode_3p	nUniqueBarcodes_5p3p	KLdivergence_5p3p	topBarcode_5p3p	secondBarcode_5p3p	thirdBarcode_5p3p	PWM_A	PWM_C	PWM_G	PWM_T
			bits = line.split("\t");
			Insert tmp = new Insert(bits[0], Integer.valueOf(bits[1]).intValue());
			if(index_5p != -1)
				tmp.addRandomBarcode5p(Integer.valueOf(bits[index_5p]).intValue(), Double.valueOf(bits[index_5p+1]).doubleValue());
			if(index_3p != -1)
				tmp.addRandomBarcode3p(Integer.valueOf(bits[index_3p]).intValue(), Double.valueOf(bits[index_3p+1]).doubleValue());
			if(index_5p3p != -1)
				tmp.addRandomBarcode5p3p(Integer.valueOf(bits[index_5p3p]).intValue(), Double.valueOf(bits[index_5p3p+1]).doubleValue());
			_allInserts.put(bits[0], tmp);
		}
		barcodeStatsReader.close();
	}*/


	/*
	 * 
	 */
	public void tidyUp() throws IOException{
		_optimalAlignmentWriter.flush();
		_optimalAlignmentWriter.close();
	}



	/**
	 * Processes the current read:
	 *  1- parses the alignments 
	 *  2- prioritises the alignment based on miRNA>tRNA>piRNA>gencode>circularRNA
	 *  3- writes accepted alignment(s)
	 *  4- [optional] looks at random barcode over-representation to moderate the read count TODO!
	 *  5- outputs the expression estimates
	 * @param thisRead
	 */
	public void assignRead(HashMap<SAMRecord, String> thisRead){
		Iterator<Entry<SAMRecord,String>> it = thisRead.entrySet().iterator();
		Entry<SAMRecord,String> tmp;

		// store the alignment records and the libraries to which they map
		HashMap<String, ArrayList<SAMRecord>> readsByLibrary = new HashMap<String, ArrayList<SAMRecord>>();  
		ArrayList<SAMRecord> keepAlignments = new ArrayList<SAMRecord>();
		boolean readHasSenseAlignment = false;

		/*
		 * Loop through alignments and count libraries that capture alignments
		 */
		int count = 0;
		while(it.hasNext()){
			tmp = it.next();
			boolean readNegativeStrand = tmp.getKey().getReadNegativeStrandFlag();
			if(!readNegativeStrand  &&  !readHasSenseAlignment)
				readHasSenseAlignment = true;

			// if this read is antisense, reverse complement it
			if(readNegativeStrand)
				tmp.getKey().setReadString(Sequence.reverseComplement(tmp.getKey().getReadString()));

			// add this read to the list of Insert sequences (only if it is the first in the alignment list)
			if(count == 0){
				String insertSeq = tmp.getKey().getReadString();
				String readID = tmp.getKey().getReadName();

				// Only add to the insert count if these have not already been counted by the random barcode summary!
				if(!_hasBarcodeStats){
					if(!_allInserts.containsKey(insertSeq))
						_allInserts.put(insertSeq, new Insert(insertSeq, 0));

					_allInserts.get(insertSeq).addRead(readID);
				}

				count++;
			}




			if(tmp.getValue().equals("miRNA")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "miRNA_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "miRNA_antisense", tmp.getKey());
			}else if(tmp.getValue().equals("tRNA")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "tRNA_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "tRNA_antisense", tmp.getKey());
			}else if(tmp.getValue().equals("piRNA")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "piRNA_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "piRNA_antisense", tmp.getKey());
			}else if(tmp.getValue().equals("gencode")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "gencode_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "gencode_antisense", tmp.getKey());
			}else if(tmp.getValue().equals("circRNA")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "circRNA_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "circRNA_antisense", tmp.getKey());
			}

			//System.out.println(tmp.getKey().getReadName()+"\ttmp.getValue()="+tmp.getValue()+"\t"+tmp.getKey().getReferenceName()+"\t"+tmp.getValue()+"\t"+(!readNegativeStrand)+"\t"+readHasSenseAlignment);
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
		//System.out.println();


		/*
		 * If this is a miRNA aligned read, check to see if it overlaps any known mature sequences
		 */
		Iterator<SAMRecord> it2 = keepAlignments.iterator();
		SAMRecord tmp2;
		if(keptLibrary.startsWith("miRNA_")){
			while(it2.hasNext()){
				tmp2 = it2.next();
				String[] refIDbits = tmp2.getReferenceName().split(":");
				String mapsTo = refIDbits[2];
				for(int i=3;i<refIDbits.length;i++)
					mapsTo = mapsTo.concat(":"+refIDbits[i]);

				// if this read maps to a mature miRNA, 
				if(miRNAprecursors.get(mapsTo).getMatureOverlapForRead(tmp2.getReadName(), tmp2.getAlignmentStart(), tmp2.getAlignmentEnd(), keptLibrary.endsWith("_antisense")).length() > 0){
					if(keptLibrary.endsWith("_sense"))
						keptLibrary = "miRNAmature_sense";
					else
						keptLibrary = "miRNAmature_antisense";
					break;
				}
			}

			// if no mature sequence explains this read, explicitly assign it to the precursor library 
			if(keptLibrary.startsWith("miRNA_")){
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
		//int counter = 0;
		while(it2.hasNext()){

			tmp2 = it2.next();
			//System.out.println(tmp2.getReferenceName());
			String insertSequence = tmp2.getReadString();

			/*
			 * Parse the SAM reference string:
			 *   [genome|nogenome]:[library]:[referenceID]
			 * e.g.
			 *   genome:gencode:ENST00000408108.1:miRNA:MIR486-201:GN=ENSG00000221035.1 
			 */
			String[] refIDbits = tmp2.getReferenceName().split(":");
			String isGenomeMapped = refIDbits[0];
			String mapsTo = refIDbits[2];
			for(int i=3;i<refIDbits.length;i++)
				mapsTo = mapsTo.concat(":"+refIDbits[i]);


			/*
			 * If this is a miRNA aligned read that HAS an overlap with a mature sequence, change the library and reference it to that of the mature
			 */
			boolean writeThisAlignment = true;
			String matureID = "NA";
			if(keptLibrary.startsWith("miRNAmature_")){
				matureID = miRNAprecursors.get(mapsTo).getMatureOverlapForRead(tmp2.getReadName(), tmp2.getAlignmentStart(), tmp2.getAlignmentEnd(), keptLibrary.endsWith("_antisense"));
				// if this read maps to a mature miRNA:
				if(matureID.length() > 0){
					_allInserts.get(insertSequence).setMapsToLib(keptLibrary);
					_allInserts.get(insertSequence).addReferenceAlignment(matureID);
				}else{
					// if read does not hit annotated mature miRNA, suppress this alignment
					writeThisAlignment = false;
				}
			}else{
				//System.out.println(insertSequence+"\t"+keptLibrary+"\t"+mapsTo);
				try{
					_allInserts.get(insertSequence).setMapsToLib(keptLibrary);
					_allInserts.get(insertSequence).addReferenceAlignment(mapsTo);
				}
				catch(NullPointerException e){
					//Thunder.printLineErr("ERROR: "+tmp2.getReadName()+"\t"+insertSequence+"\t"+keptLibrary+"\t"+mapsTo);
					ExceRpt_Tools.printLineErr("ERROR: Ignoring alignment: "+tmp2.getReadName()+"\t"+keptLibrary+"\t"+mapsTo);
					//e.printStackTrace();
				}
				//if(counter > 0  &&  !lastLib.equals(keptLibrary))
				//	System.out.println(insertSequence+"\t"+lastLib+"\t"+keptLibrary);
				//lastLib = keptLibrary;
				//counter++;
			}


			/*
			 * Write chosen alignments
			 */
			//System.out.print(tmp2.getReadName()+"\t"+mapsTo+"\t"+keptLibrary+"\t"+tmp2.getAlignmentStart()+"\t"+tmp2.getAlignmentEnd()+"\tNM:i:"+tmp2.getIntegerAttribute("NM")+"\tMD:Z:"+tmp2.getStringAttribute("MD"));
			try {
				if(writeThisAlignment)
					_optimalAlignmentWriter.write(tmp2.getReadName()+"\t"+keptLibrary+"\t"+isGenomeMapped+"\t"+mapsTo+"\t"+tmp2.getAlignmentStart()+"\t"+tmp2.getAlignmentEnd()+"\tNM:i:"+tmp2.getIntegerAttribute("NM")+"\tMD:Z:"+tmp2.getStringAttribute("MD")+"\t"+matureID+"\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}



	// Hash that holds all libraries, references, and expression statistics
	// library, referenceID, {uniqueInserts, readCount, readCount(multimap adjusted), readCount(multimap+barcode adjusted}
	HashMap<String, HashMap<String, double[]>> _library2referenceID2counts = new HashMap<String, HashMap<String, double[]>>();


	/**
	 * Loop through all inserts, create hash maps to contain expressions
	 */
	private void computeExpression(){
		//System.out.println("_allInserts.size(): "+_allInserts.size());
		Iterator<String> it = _allInserts.keySet().iterator();
		while(it.hasNext()){
			Insert tmp = _allInserts.get(it.next());
			//System.out.println("tmp.getInsertSequence(): "+tmp.getInsertSequence()+"\t"+tmp.mapsToLib());
			if(!_library2referenceID2counts.containsKey(tmp.mapsToLib()))
				_library2referenceID2counts.put(tmp.mapsToLib(), new HashMap<String, double[]>());

			int nAlignments = tmp.getReferenceAlignments().size();

			//System.out.println(tmp.getInsertSequence()+"\t"+tmp.getInsertSequenceCount()+"\t"+tmp.mapsToLib()+"\t"+nAlignments);

			for(String thisAlignment: tmp.getReferenceAlignments()){
				double[] alignmentCounts;
				if(!_library2referenceID2counts.get(tmp.mapsToLib()).containsKey(thisAlignment))
					alignmentCounts = new double[]{0.0,0.0,0.0,0.0};
				else
					alignmentCounts = _library2referenceID2counts.get(tmp.mapsToLib()).get(thisAlignment);

				alignmentCounts[0] += 1;
				alignmentCounts[1] += tmp.getInsertSequenceCount();
				alignmentCounts[2] += (tmp.getInsertSequenceCount()+0.0)/(nAlignments+0.0);
				alignmentCounts[3] += (tmp.getBarcodeCount()+0.0)/(nAlignments+0.0);

				_library2referenceID2counts.get(tmp.mapsToLib()).put(thisAlignment, alignmentCounts);
			}
		}
		//System.out.println("_library2referenceID2counts.size() = "+_library2referenceID2counts.size());
	}


	/**
	 * Write the insert counts, read counts, and adjusted counts for each of the observed libraries 
	 * @param basePath
	 * @throws IOException
	 */
	public void writeCounts(String basePath) throws IOException{
		Iterator<String> libraryIterator = _library2referenceID2counts.keySet().iterator();
		while(libraryIterator.hasNext()){
			String thisLibrary = libraryIterator.next();
			//System.out.println(thisLibrary+"\t"+_library2referenceID2counts.get(thisLibrary).size());

			if(_library2referenceID2counts.get(thisLibrary).size() > 0){
				BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/readCounts_"+thisLibrary+".txt"));
				out.write("ReferenceID\tuniqueReadCount\ttotalReadCount\tmultimapAdjustedReadCount\tmultimapAdjustedBarcodeCount\n");

				TreeMap<String,double[]> tmp_sorted;
				tmp_sorted = new TreeMap<String,double[]>(new Comparator_SortByAdjustedReadCount(_library2referenceID2counts.get(thisLibrary)));
				tmp_sorted.putAll(_library2referenceID2counts.get(thisLibrary));
				Iterator<String> it = tmp_sorted.keySet().iterator();
				while(it.hasNext()){
					String thisID = it.next();
					double[] quants = _library2referenceID2counts.get(thisLibrary).get(thisID);
					out.write(thisID+"\t"+Math.round(quants[0])+"\t"+Math.round(quants[1])+"\t"+quants[2]+"\t"+quants[3]+"\n");
				}
				out.flush();
				out.close();


				// If this is the gencode library, compute gene level expressions (as this is slow to do in post-processing 
				if(thisLibrary.startsWith("gencode_")){

					out = new BufferedWriter(new FileWriter(basePath+"/readCounts_"+thisLibrary+"_geneLevel.txt"));
					out.write("ReferenceID\tuniqueReadCount\ttotalReadCount\tmultimapAdjustedReadCount\tmultimapAdjustedBarcodeCount\n");

					HashMap<String, double[]> geneLevelQuants = new HashMap<String, double[]>(); 
					it = _library2referenceID2counts.get(thisLibrary).keySet().iterator();
					while(it.hasNext()){
						String thisID = it.next();
						//System.out.println("thisID = "+thisID);
						/*   [library]:[referenceID]
						 *   ENST00000431311.1:lincRNA:RP11-214L19.1-001:GN=ENSG00000237290.1  */
						String[] refIDbits = thisID.split(":");
						if(refIDbits.length < 4){
							// TODO: figure this out!!!
							//System.err.println("ERROR parsing "+thisID+" from library "+thisLibrary);
						}else{
							String[] geneNameBits = refIDbits[2].split("-");
							String geneName = geneNameBits[0];
							for(int i=1;i<geneNameBits.length-1;i++)
								geneName.concat("-"+geneNameBits[i]);
							//String newID = geneName+":"+refIDbits[3].substring(3)+":"+refIDbits[1];
							String newID = geneName+":"+refIDbits[1];

							if(!geneLevelQuants.containsKey(newID))
								geneLevelQuants.put(newID, new double[]{0.0,0.0,0.0,0.0});
							double[] thisGeneCounts = geneLevelQuants.get(newID);

							double[] thisTranscriptCounts = _library2referenceID2counts.get(thisLibrary).get(thisID);
							if(thisTranscriptCounts[0] > thisGeneCounts[0])
								thisGeneCounts[0] = thisTranscriptCounts[0];
							if(thisTranscriptCounts[1] > thisGeneCounts[1])
								thisGeneCounts[1] = thisTranscriptCounts[1];
							thisGeneCounts[2] += thisTranscriptCounts[2];
							thisGeneCounts[3] += thisTranscriptCounts[3];

							geneLevelQuants.put(newID, thisGeneCounts);
						}
					}

					tmp_sorted = new TreeMap<String,double[]>(new Comparator_SortByAdjustedReadCount(geneLevelQuants));
					tmp_sorted.putAll(geneLevelQuants);
					it = tmp_sorted.keySet().iterator();
					while(it.hasNext()){
						String thisID = it.next();
						double[] quants = geneLevelQuants.get(thisID);
						out.write(thisID+"\t"+Math.round(quants[0])+"\t"+Math.round(quants[1])+"\t"+quants[2]+"\t"+quants[3]+"\n");
					}
					out.flush();
					out.close();
				}
			}
		}
	}



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



	/**
	 * Reads the read alignments
	 * reference ID format is
	 *   [genome|nogenome]:[library]:[referenceID]
	 * e.g.
	 *   genome:miRNA:hsa-mir-486-1:MI0002470:Homo:sapiens:miR-486:stem-loop
	 * 
	 * @param path_readAlignments
	 */
	public boolean read_Reads(File path_readAlignments){
		/*
		 * Read the read alignments
		 */
		SAMFileReader inputSam = new SAMFileReader(path_readAlignments);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		//inputSam.setValidationStringency(ValidationStringency.LENIENT);

		// TODO: put this if/else back when we're done testing
		//if(inputSam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.queryname)){
		SAMRecord thisRecord;
		HashMap<SAMRecord, String> thisRead = new HashMap<SAMRecord, String>(); 
		SAMRecordIterator it = inputSam.iterator();
		String lastReadID = null;
		while(it.hasNext()){
			thisRecord = it.next();
			//System.out.println(thisRecord.getReferenceName());
			if(!thisRecord.getReadName().equals(lastReadID)  &&  lastReadID != null){
				// new, non first
				assignRead(thisRead);
				thisRead = new HashMap<SAMRecord, String>();
			}

			// put the SAM record into the map with the library type as the value
			thisRead.put(thisRecord, thisRecord.getReferenceName().split(":")[1]);
			lastReadID = thisRecord.getReadName();
			//System.out.println(thisRecord.getReferenceName().split(":")[0]+"\t"+thisRecord.getReferenceName()+"\t"+thisRecord.getReadName());

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
		options.addOption(OptionBuilder.withArgName(".SAM or .BAM").hasArg().withDescription("Path to READ alignments to the hairpins").create("reads2all"));
		options.addOption(OptionBuilder.withArgName(".stats").hasArg().withDescription("[optional] Path to random barcode stats").create("randombarcode"));
		options.addOption(OptionBuilder.withArgName("csv list").hasArg().withDescription("[optional] Library priorities for quantification. Comma separated list of libraries in *descending* order of importance. Default: miRNA,tRNA,piRNA,gencode,circRNA - this can also be used to suppress libraries during quantification.").create("libPriority"));
		options.addOption(OptionBuilder.withArgName("directory").hasArg().withDescription("Base path to write the results into").create("outputPath"));
		return options;
	}


	public static void main(String[] args) throws ParseException, IOException {
		/*String hairpin2genome = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MICRO_RNA/miRBase_v21_hairpin_hsa_hg19_aligned.sam";
		String mature2hairpin = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MICRO_RNA/miRBase_v21_mature_hairpin_hsa_aligned.sam";
		//String reads_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_ALL.sam";
		String reads_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_LIBS.sam";
		//String reads_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_top.sam";
		//String reads_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/test2.sam";
		String output_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants";
		args = new String[]{"ProcessEndogenousAlignments",
				"--hairpin2genome",hairpin2genome,
				"--mature2hairpin",mature2hairpin,
				"--reads2all",reads_path,
				"--outputPath",output_path,
				"--libPriority","miRNA,tRNA,piRNA,gencode,circRNA"
		};*/


		CommandLine cmdArgs = ExceRpt_Tools.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption("hairpin2genome") && cmdArgs.hasOption("mature2hairpin") && cmdArgs.hasOption("reads2all") && cmdArgs.hasOption("outputPath")){
			//ProcessEndogenousAlignments engine = new ProcessEndogenousAlignments(new File(cmdArgs.getOptionValue("hairpin2genome")), new File(cmdArgs.getOptionValue("mature2hairpin")));
			ProcessEndogenousAlignments engine = new ProcessEndogenousAlignments(cmdArgs.getOptionValue("outputPath"));

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

			// Read hairpin alignments to the genome and mature alignments to the hairpins
			ExceRpt_Tools.printLineErr("Reading miRNA annotation info");
			engine.read_miRNAinfo(new File(cmdArgs.getOptionValue("hairpin2genome")), new File(cmdArgs.getOptionValue("mature2hairpin")));

			// read and process the RNA-seq read alignments
			ExceRpt_Tools.printLineErr("Processing RNA-seq alignments");
			boolean cont = engine.read_Reads(new File(cmdArgs.getOptionValue("reads2all")));

			if(cont){
				ExceRpt_Tools.printLineErr("Computing expressions");
				engine.computeExpression();

				// write the miRNA counts
				ExceRpt_Tools.printLineErr("Writing read counts and tidying up");
				//engine.writeCounts_miRNA(cmdArgs.getOptionValue("outputPath"));
				// write the counts from other libraries
				engine.writeCounts(cmdArgs.getOptionValue("outputPath"));

				// close global buffered writer(s)
				engine.tidyUp();

				ExceRpt_Tools.printLineErr("Done!");
			}
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(ExceRpt_Tools.ExceRpt_Tools_EXE_COMMAND+" ProcessEndogenousAlignments", getCmdLineOptions());
			System.err.println();
		}
	}
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


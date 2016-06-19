package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class QuantifyEndogenousAlignments {

	private BufferedReader _alignmentReader;
	public QuantifyEndogenousAlignments(String inputPath) throws FileNotFoundException{
		_alignmentReader = new BufferedReader(new FileReader(inputPath));
	}


	private boolean _hasBarcodeStats = false;
	//private HashMap<String, Insert> _allInserts = new HashMap<String, Insert>();



	// Hash that holds all libraries, references, and expression statistics
	// library, referenceID, {uniqueInserts, readCount, readCount(multimap adjusted), readCount(multimap+barcode adjusted}
	HashMap<String, HashMap<String, double[]>> _library2referenceID2counts = new HashMap<String, HashMap<String, double[]>>();





	/**
	 * Loop through all inserts, create hash maps to contain expressions
	 */
	private void computeExpression(){
		//System.out.println("_allInserts.size(): "+_allInserts.size());
		/*Iterator<String> it = _allInserts.keySet().iterator();
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
		}*/
		//System.out.println("_library2referenceID2counts.size() = "+_library2referenceID2counts.size());
	}





	/**
	 * Write the insert counts, read counts, and adjusted counts for each of the observed libraries 
	 * @param basePath
	 * @throws Exception 
	 */
	public void writeCounts(String basePath, String annotationPath) throws Exception{
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

					//read gencode gtf and collapse to gene id
					//TranscriptAnnotation annotation = ReadGTF.readGTF(annotationPath);
					
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
						if(refIDbits.length < 3){
							// TODO: figure this out!!!
							//System.err.println("ERROR parsing "+thisID+" from library "+thisLibrary);
						}else{
							//String geneID = annotation.getGeneForTranscript(refIDbits[0]);
							//annotation.getTranscript(refIDbits[0]).get
							
							String[] geneNameBits = refIDbits[2].split("-");
							String geneName = geneNameBits[0];
							for(int i=1;i<geneNameBits.length-1;i++)
								geneName = geneName.concat("-"+geneNameBits[i]);
							//String newID = geneName+":"+refIDbits[3].substring(3)+":"+refIDbits[1];
							//String newID = geneName+":"+refIDbits[1];
							String newID = geneName;

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

	private HashMap<Integer, String> _alignmentDictionary = new HashMap<Integer, String>();
	/**
	 * Reads the map between alignment index and alignment string (to save space and memory where there are a lot of different alignments)
	 * @param path
	 * @throws IOException
	 */
	private void readAlignmentDictionary(String path) throws IOException{
		BufferedReader reader = new BufferedReader(new FileReader(new File(path)));
		String line;
		String[] bits;
		while((line=reader.readLine())!=null){
			bits = line.trim().split("\t");
			_alignmentDictionary.put(Integer.valueOf(bits[0]), bits[1]);
		}
		reader.close();
	}
	

	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName(".txt").hasArg().withDescription("Path to accepted READ alignments output by 'ProcessEndogenousAlignments'").create("acceptedAlignments"));
		options.addOption(OptionBuilder.withArgName(".dict").hasArg().withDescription("Path to the dictionary of alignments output by 'ProcessEndogenousAlignments'").create("dict"));
		//options.addOption(OptionBuilder.withArgName(".stats").hasArg().withDescription("[optional] Path to random barcode stats").create("randombarcode"));
		options.addOption(OptionBuilder.withArgName("directory").hasArg().withDescription("Base path to write the results into").create("outputPath"));
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("Path to the GTF file containing the transcript/gene relationship").create(Thunder.OPT_PATH_ANNOTATION));
		return options;
	}



	public static void main(String[] args) throws Exception {
		/*
		//String input_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_Accepted.sorted.txt";
		String input_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/acceptedAlignments_NEW.sorted.txt";		
		//String input_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_Accepted.sorted.head.txt";
		//String input_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_Accepted.sorted.head10.txt";
		String output_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/new";
		args = new String[]{"QuantifyEndogenousAlignments",
				"--acceptedAlignments",input_path,
				"--outputPath",output_path,
		};*/

		CommandLine cmdArgs = ExceRpt_Tools.parseArgs(args, getCmdLineOptions());

		//if(cmdArgs.hasOption("acceptedAlignments") && cmdArgs.hasOption("outputPath")){
		if(cmdArgs.hasOption("acceptedAlignments") && cmdArgs.hasOption("outputPath") && cmdArgs.hasOption("dict")){// && cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION)){

			ExceRpt_Tools.printLineErr("Reading alignments in: "+cmdArgs.getOptionValue("acceptedAlignments"));
			ExceRpt_Tools.printLineErr("Using alignment dictionary: "+cmdArgs.getOptionValue("dict"));
			ExceRpt_Tools.printLineErr("Writing results to: "+cmdArgs.getOptionValue("outputPath"));

			try {
				QuantifyEndogenousAlignments engine = new QuantifyEndogenousAlignments(cmdArgs.getOptionValue("acceptedAlignments"));

				ExceRpt_Tools.printLineErr("Reading alignment dictionary...");
				engine.readAlignmentDictionary(cmdArgs.getOptionValue("dict"));
				
				ExceRpt_Tools.printLineErr("Reading alignments...");
				engine.readAndCountInserts();

				ExceRpt_Tools.printLineErr("Writing read counts...");
				engine.writeCounts(cmdArgs.getOptionValue("outputPath"), cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION));
				
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}



			ExceRpt_Tools.printLineErr("Done!");

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(ExceRpt_Tools.ExceRpt_Tools_EXE_COMMAND+" QuantifyEndogenousAlignments", getCmdLineOptions());
			System.err.println();
		}
	}




	/**
	 * Read all reads and alignments corresponding to a given insert sequence
	 * @param thisInsert
	 */
	private void addInsert(ArrayList<String> thisInsert){
		Iterator<String> it;
		String[] tmpBits = thisInsert.get(0).trim().split("\t");

		String insertSequence = tmpBits[1];
		String library = tmpBits[2];
		
		HashSet<String> referenceIDs = new HashSet<String>();
		String[] alignmentBits = _alignmentDictionary.get(Integer.valueOf(tmpBits[4])).split("\\|");
		for(String x:alignmentBits)
			referenceIDs.add(x.trim());
		
		HashSet<String> readIDs = new HashSet<String>();
		

		// Get the unique set of read IDs and reference alignments for this insert
		it = thisInsert.iterator();
		while(it.hasNext()){
			String[] thisAlignment = it.next().trim().split("\t");
			//System.out.println("--"+thisAlignment[0]+"\t"+thisAlignment[1]+"\t"+thisAlignment[2]+"\t"+thisAlignment[4]);
			readIDs.add(thisAlignment[0]);
			/*if(library.equals("miRNAmature_sense")  ||  library.equals("miRNAmature_antisense"))
				referenceIDs.add(thisAlignment[9]);
			else
				referenceIDs.add(thisAlignment[4]);*/

			/*String[] alignmentBits = thisAlignment[4].split("\\|");
			for(String x:alignmentBits){
				//System.out.println(x);
				referenceIDs.add(x.trim());
			}*/
			//System.out.println(alignmentBits.length);

		}

		Insert insert = new Insert(insertSequence, 0);
		it = readIDs.iterator();
		while(it.hasNext())
			insert.addRead(it.next());


		for(String thisReferenceID: referenceIDs){
			double[] alignmentCounts;
			if(!_library2referenceID2counts.containsKey(library))
				_library2referenceID2counts.put(library, new HashMap<String, double[]>());

			if(!_library2referenceID2counts.get(library).containsKey(thisReferenceID))
				alignmentCounts = new double[]{0.0,0.0,0.0,0.0};
			else
				alignmentCounts = _library2referenceID2counts.get(library).get(thisReferenceID);

			alignmentCounts[0] += 1;
			alignmentCounts[1] += insert.getInsertSequenceCount();
			alignmentCounts[2] += (insert.getInsertSequenceCount()+0.0)/(referenceIDs.size()+0.0);
			alignmentCounts[3] += (insert.getBarcodeCount()+0.0)/(referenceIDs.size()+0.0);

			_library2referenceID2counts.get(library).put(thisReferenceID, alignmentCounts);
		}

		//System.out.println("nInserts: "+insert.getInsertSequenceCount()+"\tnBarcodes: "+insert.getBarcodeCount());
	}



	/**
	 * 
	 * @throws IOException
	 */
	private void readAndCountInserts() throws IOException{
		ArrayList<String> thisInsert = new ArrayList<String>();
		int count = 0;

		try{
			while((_thisLine=_alignmentReader.readLine()) != null){
				//System.out.println(">"+_thisLine);
				_thisInsert = _thisLine.trim().split("\t")[1]; 
				if(_lastInsert == null  ||  _thisInsert.equals(_lastInsert)){ // this is the same insert sequence, add it and keep going!
					thisInsert.add(_thisLine);
					_lastInsert = _thisInsert;
				}else{ // this is a new insert, keep the line, but return the current ArrayList
					addInsert(thisInsert);
					count ++;
					thisInsert = new ArrayList<String>();

					thisInsert.add(_thisLine);
					_lastInsert = _thisInsert;
				}
			}
			// add final insert:
			addInsert(thisInsert);
		}catch(Exception e){
			System.err.println("ERROR: failed on line:\n"+_thisLine+"\n\n");
			e.printStackTrace();
		}

		/*while((thisInsert=readNextInsertSeq()).size() > 0){
			System.out.println("\n"+count+"\t"+thisInsert.size());
			addInsert(thisInsert);
			count ++;
			if(_doneReading)
				break;
		}*/
		ExceRpt_Tools.printLineErr(" - Read alignments for "+(count+1)+" insert sequences");
	}


	private String _thisLine = "";
	private String _lastInsert = null, _thisInsert = "";

	//private boolean _doneReading = false;
	/*public ArrayList<String> readNextInsertSeq() throws IOException{
		ArrayList<String> newInsert = new ArrayList<String>();
		if(_lastLine != null) // this is NOT the first time this subroutine has been called
			newInsert.add(_lastLine);
			return(readNextInsertSeq(newInsert));
	}

	public ArrayList<String> readNextInsertSeq(ArrayList<String> thisInsert) throws IOException{
		if((_thisLine=_alignmentReader.readLine()) != null){
			System.out.println(">"+_thisLine);
			_thisInsert = _thisLine.trim().split("\t")[1]; 
			if(_lastInsert == null  ||  _thisInsert.equals(_lastInsert)){ // this is the same insert sequence, add it and keep going!
				thisInsert.add(_thisLine);
				_lastLine = _thisLine;
				_lastInsert = _thisInsert;
				return readNextInsertSeq(thisInsert);
			}else{ // this is a new insert, keep the line, but return the current ArrayList
				_lastLine = _thisLine;
				_lastInsert = _thisInsert;
				return thisInsert;
			}
		}else{ // there are no more lines, just return the current list of inserts
			_doneReading = true;
			return thisInsert;
		}
	}
	 */


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
}







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




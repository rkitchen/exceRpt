package main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import utils.IO_utils;
import annotation.NCBITaxonomy_Engine;
import annotation.NCBITaxonomy_Node;

public class ProcessExogenousAlignments {

	private NCBITaxonomy_Engine _taxonomy;
	private boolean _topDownSearch = false;
	private double _minFraction = 0.0;

	/**
	 * Main
	 * 
	 * @param taxonomyPath
	 * @param alignmentsPath
	 * @throws IOException
	 */
	public ProcessExogenousAlignments(int dataType, boolean topDown, double minFraction, int minReads, int verbose) throws IOException{ 
		_verbose = verbose; 
		_minReads = minReads;
		_topDownSearch = topDown;
		_minFraction = minFraction;
		_thisDataType = dataType;

		//System.err.println("_verbose = "+_verbose+"\t_verbose == 2: "+(_verbose == 2));

		// if these alignments are from RDP, we have to to top-down as the leaf nodes are not guaranteed to be the ones with the read alignments
		//if(dataType == DATA_TYPE_RRNA)
		//	_topDownSearch = true;
	}


	/**
	 * 
	 * @param taxonomyPath
	 * @throws IOException
	 */
	public void readTaxonomy(String taxonomyPath) throws IOException{
		_taxonomy = new NCBITaxonomy_Engine(taxonomyPath);
	}



	private BufferedReader _alignmentReader;

	/**
	 * 
	 * @param alignmentsPath
	 * @throws IOException
	 */
	public void processAlignments(String alignmentsPath, int batchSize) throws IOException{
		//IO_utils.printLineErr("Reading alignments...");
		_alignmentReader = new BufferedReader(new FileReader(alignmentsPath));
		int batchCount = 1;
		//int readCounter = 0;
		//boolean hasMoreAlignments = false;
		//System.out.println("batch "+batchCount);
		while(readAlignments(batchSize)){
			// printing:
			/*System.out.println("batch "+batchCount);
			Iterator<String> it = _reads2species.keySet().iterator();
			while(it.hasNext())
				System.out.println("\t"+it.next());
			 */

			int inputReadNumber = _reads2species.size();
			int minReads = (int)Math.floor(_minFraction * (inputReadNumber+0.0));


			IO_utils.printLineErr("Processing alignments in batch "+batchCount);
			//System.err.println(" (batch "+batchCount+")");
			if(_countIgnoredReads)
				IO_utils.printLineErr(" - Ignored "+_ignoredReads.size()+" alignments to species that could not be found in the taxonomy ("+countIgnoredReads()+" reads removed entirely)");
			IO_utils.printLineErr(" - Detected "+_species2reads.size()+" possible species from "+inputReadNumber+" aligned reads (minReads="+minReads+")");
			//readCounter += _reads2species.size();




			if(_topDownSearch){
				// top down - slow but complete
				//processAlignments(_taxonomy.getRootNode());
				processAlignments(_taxonomy.getNode("cellular organisms"));
				//processAlignments(_taxonomy.getNode("bacteria"));
			}else{
				// bottom up - fast (default)
				processAlignments_bottomUP(new HashSet<String>(), minReads);
				IO_utils.printLineErr(" - Done. Suppressed "+_reads2species.size()+" total reads ("+(Math.round((_reads2species.size()+0.0)*100000.0/(inputReadNumber+0.0))/1000.0)+"% of input reads) due to minimum read requirement.");
			}


			//IO_utils.printLineErr(" - Assigned "+_readAssignment.size()+" of "+readCounter+" reads ("+Math.round((_readAssignment.size()+0.0)*1000.0/(readCounter+0.0))/10.0+"%) to a subset of the taxonomy");
			batchCount ++;

			// clean up ready for the next batch
			_reads2species = new HashMap<String, HashSet<String>>();
			_species2reads = new HashMap<String, HashSet<String>>();

			//if(hasMoreAlignments)
			//	IO_utils.printLineErr("Reading more alignments...");
		}

		IO_utils.printLineErr("Quantifying");
		countReads();
		//IO_utils.printLineErr("Quantified taxonomy:");
		String indent = ">";
		_taxonomy.getRootNode().printSummary(indent, _minReads);

		IO_utils.printLineErr("Done");
	}


	private String _lastReadID=null, _lastSpeciesID=null;

	private int _minReads = 3;

	private int _thisDataType = 0;
	private static final int DATA_TYPE_GENOME = 0;
	private static final int DATA_TYPE_RRNA = 1;

	/**
	 * 
	 * @param alignmentsPath where to 
	 * @param batchSize number of reads to read
	 * 
	 * @return true if there are more alignments to read
	 * @throws IOException
	 */
	private boolean readAlignments(int batchSize) throws IOException{

		IO_utils.printLineErr("Reading alignments...");

		boolean hasMoreAlignments = false;
		int readCounter = 0;
		String line;

		if(_lastReadID != null){// there is a read left over from the last batch, add it!
			addRead(_lastReadID, _lastSpeciesID);
		}

		//int count = 0;
		while((line=_alignmentReader.readLine())!=null){
			hasMoreAlignments = true;

			String readID = "";
			String species = "";

			if(_thisDataType == DATA_TYPE_GENOME){
				String[] bits = line.split("\t");
				readID = bits[0].trim();
				//String species = bits[2].replace("_", " ").toLowerCase().trim();
				species = NCBITaxonomy_Engine.sanitiseName(bits[2]);				
			}
			else if(_thisDataType == DATA_TYPE_RRNA){
				String[] bits = line.split(" ");
				readID = bits[0].trim();
				species = findMatchToTaxonomy(bits[1].trim());
			}

			//if(count < 10){
			//System.out.println(_lastReadID+"\t"+readID+"\t"+line);
			//}
			if(readID.equals(_lastReadID)  ||  _lastReadID == null){
				addRead(readID, species);
				_lastReadID = readID;
				_lastSpeciesID = species;
			}else{ // new read
				readCounter ++;
				//System.out.println("readCounter: "+readCounter);
				_lastReadID = readID;
				_lastSpeciesID = species;
				if(readCounter < batchSize)
					addRead(readID, species);
				else
					break;
			}
		}

		if(!hasMoreAlignments)
			_lastReadID = null;



		//System.out.println(_maxCountSpecies+": "+_maxCount);
		//br.close();

		return hasMoreAlignments;
	}



	/**
	 * Specialist function to search for the appropriate NCBI taxonomy ID from the provided RDP ID.
	 * 
	 * @param speciesStringRDP
	 * @return
	 */
	private String findMatchToTaxonomy(String speciesStringRDP){
		//boolean matchedID = false;
		String matchedTaxaID = null;

		// split the RDP ID into the main and lineage parts
		String[] mainBits = speciesStringRDP.split("Lineage=");
		String firstPart = NCBITaxonomy_Engine.sanitiseName(mainBits[0]);
		//String firstPart = mainBits[0].replaceAll("_", " ").replaceAll("\"", "").toLowerCase().trim();
		firstPart = firstPart.split(";")[0].trim();
		//System.out.println("firstPart = "+firstPart);

		// look for a taxa match to the first part of the RDP ID (preferred)
		if(_taxonomy.containsNode(firstPart)){
			// if the main ID is in the taxonomy, use it:
			//System.out.println(firstPart+"\t"+_taxonomy.containsNode(firstPart));
			matchedTaxaID = firstPart;
		}

		// if the main ID is NOT in the taxonomy, try different substrings of it, e.g. remove the first section that is usually an RDP specific ID:
		if(matchedTaxaID == null){
			String[] tmpBits2 = firstPart.split(" ");
			if(tmpBits2.length >= 2){
				// try just removing the first part of the ID
				String toTry = tmpBits2[1];
				for(int i=2;i<tmpBits2.length;i++)
					toTry = toTry.concat(" "+tmpBits2[i]);
				//System.out.println(">>> "+toTry+"\t"+_taxonomy.containsNode(toTry));
				if(_taxonomy.containsNode(toTry)){
					matchedTaxaID = toTry;
				}
			}
		}


		if(matchedTaxaID != null){
			//System.out.println("FOUND:\t"+matchedTaxaID+"\t"+_taxonomy.getNode(matchedTaxaID).getLevel());
		}


		//
		// FINALLY just check the Lineage string in case there is a better description here (search string from end to start to find the best possible match)
		//
		String matchedID_lineage = null;
		String[] tmpBits = mainBits[1].split(";");
		for(int i=tmpBits.length-1;i>0;i--){
			String tmpID = tmpBits[i].toLowerCase();
			tmpID = tmpID.replaceAll("\"", "");
			tmpID = tmpID.replaceAll("_", " ");
			if(_taxonomy.containsNode(tmpID)){
				//System.out.println(tmpID+"\t"+_taxonomy.containsNode(tmpID));
				matchedID_lineage = tmpID;
				break;
			}
		}

		if(matchedID_lineage != null){
			//System.out.println("FOUND:\t"+matchedID_lineage+"\t"+_taxonomy.getNode(matchedID_lineage).getLevel());
		}


		// If the bacteria is just 'uncultured_bacterium', pick the lineage version instead
		if(matchedTaxaID == null  ||  matchedTaxaID.equals("uncultured bacterium"))
			matchedTaxaID = matchedID_lineage;


		// TODO: compare matchedID_lineage with matchedTaxaID and return the one that is most specific

		//if(matchedTaxaID != null)
		return matchedTaxaID;
		//else
		//	return matchedID_lineage;
	}



	private boolean _countIgnoredReads = false; // can cause high memory usage if there are a lot of these
	HashSet<String> _ignoredReads = new HashSet<String>();
	/**
	 * 
	 * @param readID
	 * @param species
	 */
	private void addRead(String readID, String species){
		// check that this species is in the taxonomy
		//if(_taxonomy.getNodeName2nodeIndex().containsKey(species)){
		if(_taxonomy.containsNode(species)){

			if(_taxonomy.getNode(species).hasChildren()){
				// if we're adding a read to a non-leaf node, add it to all leaves below this node

				//System.out.println("N children = "+_taxonomy.getNode(species).getAllLeafSpeciesForNode().size());
				for(String tmp: _taxonomy.getNode(species).getAllLeafSpeciesForNode()){
					addRead(readID, tmp);
				}
			}else{
				// add the species if we haven't seen it before
				if(!_species2reads.containsKey(species)){
					_species2reads.put(species, new HashSet<String>());
					_species2readCounts.put(species, 0);
				}

				// add the read if we haven't seen it before
				if(!_reads2species.containsKey(readID)){
					_reads2species.put(readID, new HashSet<String>());
				}
				_reads2species.get(readID).add(species);

				// if this read is not already assigned to this species, add it and increment the count
				if(!_species2reads.get(species).contains(readID)){
					_species2reads.get(species).add(readID);
					int thisCount = _species2readCounts.get(species)+1;
					_species2readCounts.put(species, thisCount);
				}

			}



			//System.out.println(species+": "+_species2readCounts.get(species));
		}else if(_countIgnoredReads){
			//IO_utils.printLineErr("\'"+species+"\' not found in the taxonomy, ignoring");
			_ignoredReads.add(readID);
			//ignoredSpecies.add(species);
		}
	}
	private int countIgnoredReads(){
		int ignoredReadCount = 0;
		// count reads that are removed completely due to taxonomy ambiguity
		Iterator<String> it = _ignoredReads.iterator();
		while(it.hasNext()){
			if(!_reads2species.containsKey(it.next()))
				ignoredReadCount ++;
		}
		return ignoredReadCount;
	}






	/**
	 * recount reads so as to be cumulative going up the taxonomic tree 
	 */
	public void collapseReadCounts(){


	}


	/**
	 * Count assigned reads 
	 */
	@SuppressWarnings("unchecked")
	private void countReads(){
		Iterator<String> it = _readAssignment.keySet().iterator();
		while(it.hasNext()){
			NCBITaxonomy_Node node = _readAssignment.get(it.next());
			if(!_nodeReadCount.containsKey(node))
				_nodeReadCount.put(node, 0);
			_nodeReadCount.put(node, _nodeReadCount.get(node)+1);
			_taxonomy.getNode(node.getNodeIndex()).addRead("");
		}

		Object[] a = _nodeReadCount.entrySet().toArray();
		Arrays.sort(a, new Comparator<Object>() {
			public int compare(Object o1, Object o2) {
				return ((Map.Entry<NCBITaxonomy_Node, Integer>) o2).getValue().compareTo(
						((Map.Entry<NCBITaxonomy_Node, Integer>) o1).getValue());
			}
		});
		/*for (Object e : a) {
	    	NCBITaxonomy_Node tmp = ((Map.Entry<NCBITaxonomy_Node, Integer>) e).getKey();
	        System.out.println(((Map.Entry<NCBITaxonomy_Node, Integer>) e).getValue()+"\t"+tmp.getTaxonomyName()+"\t"+tmp.getLevel()+"");
	    }*/
	}


	private double _minReaadFrac = 0.98;
	public void setMinReadFractionForTreeDescent(double minFraction){
		_minReaadFrac = minFraction;
	}


	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void processAlignments_bottomUP(HashSet<String> assignedReadsInThisBatch, int minReads){
		// 1: sort species by # mapped reads
		Object[] a = _species2readCounts.entrySet().toArray();
		Arrays.sort(a, new Comparator() {
			public int compare(Object o1, Object o2) {
				return ((Map.Entry<String, Integer>) o2).getValue().compareTo(
						((Map.Entry<String, Integer>) o1).getValue());
			}
		});

		String mostAbundantSpecies = ((Map.Entry<String, Integer>) a[0]).getKey();

		if(_species2readCounts.get(mostAbundantSpecies) >= minReads){

			if(_verbose > 0)
				IO_utils.printLineErr(" - Reads remaining: "+_reads2species.size()+"\tResolving alignments to "+mostAbundantSpecies+" ("+_species2readCounts.get(mostAbundantSpecies)+" reads)");
			//IO_utils.printLineErr(" - Reads remaining: "+_reads2species.size()+"\tResolving alignments to "+mostAbundantSpecies+" ("+_species2reads.get(mostAbundantSpecies).size()+" reads)");

			if(_species2reads.containsKey(mostAbundantSpecies)  &&  _species2reads.get(mostAbundantSpecies).size() > 0)
				assessNode(_taxonomy.getNode(mostAbundantSpecies), _species2reads.get(mostAbundantSpecies), assignedReadsInThisBatch);
			else{
				IO_utils.printLineErr(" WARNING: failed to process reads for species:"+mostAbundantSpecies+ "(species2readCounts="+_species2readCounts.get(mostAbundantSpecies)+", in species2reads="+_species2reads.containsKey(mostAbundantSpecies));
				_species2readCounts.put(mostAbundantSpecies, 0);
			}
			//else
			//	assessNode(_taxonomy.getNode("cellular organisms"), _species2reads.get(mostAbundantSpecies), assignedReadsInThisBatch);
			//IO_utils.printLineErr("Reads remaining: "+_reads2species.size());

			if(_reads2species.size() > 0)
				processAlignments_bottomUP(assignedReadsInThisBatch, minReads);
		}

	}

	private void assessNode(NCBITaxonomy_Node thisNode, HashSet<String> reads, HashSet<String> assignedreads){
		//NCBITaxonomy_Node parent = thisNode.getParent();
		//HashSet<String> leafNodes = parent.getAllLeafSpeciesForNode();
		//System.out.println("parent: "+parent.getName()+"\tN leaf nodes: "+leafNodes.size());
		HashSet<String> leafNodes = thisNode.getAllLeafSpeciesForNode();

		if(_verbose == 2)
			IO_utils.printErr("parent: "+thisNode.getName()+"\tlevel: "+thisNode.getLevel()+"\tN leaf nodes: "+leafNodes.size());

		// loop over all reads, asking whether this parent can explain each set of alignments
		Iterator<String> it;
		//Iterator<String> reads = _reads2species.keySet().iterator();
		Iterator<String> readIterator = reads.iterator();
		while(readIterator.hasNext()){
			String thisRead = readIterator.next();

			int alignedSpeciesUnderCurrentParent = 0;
			if(_reads2species.containsKey(thisRead)){
				HashSet<String> theseSpecies = _reads2species.get(thisRead);
				it = theseSpecies.iterator();
				while(it.hasNext())
					if(leafNodes.contains(it.next()))
						alignedSpeciesUnderCurrentParent ++;

				double frac = (alignedSpeciesUnderCurrentParent+0.0)/(theseSpecies.size()+0.0);

				//System.out.println(thisRead+"\t"+thisChild.getName()+" ("+thisChild.getLevel()+"): "+frac);

				if(frac >= _minReaadFrac){ // this read works at this level
					_readAssignment.put(thisRead, thisNode);
					assignedreads.add(thisRead);
					//System.out.println(thisRead+"\t"+thisChild.getName()+" ("+thisChild.getLevel()+"): "+frac);
				}
			}
		}

		reads = tidyUp(reads, assignedreads);

		// try the parent of this node		
		//if(!thisNode.getParent().getName().equals("root"))
		//if(thisNode.getParent().getName().equals("root")  ||  reads.size() == 0){
		try{
			if(thisNode.getParent().getName().equals("cellular organisms")  ||  reads.size() == 0){
				// if we've not yet assigned these reads, assign them all to cellular organisms
				NCBITaxonomy_Node junkNode = _taxonomy.getNode("cellular organisms");

				// if we've not yet assigned these reads, assign them all to cellular organisms
				//NCBITaxonomy_Node junkNode = _taxonomy.getNode("root");

				if(_verbose == 2)
					IO_utils.printErr("parent: "+junkNode.getName()+"\tlevel: "+junkNode.getLevel());
				readIterator = reads.iterator();
				while(readIterator.hasNext()){
					String thisRead = readIterator.next();
					_readAssignment.put(thisRead, junkNode);
					assignedreads.add(thisRead);
				}
				tidyUp(reads, assignedreads);
			}else{
				assessNode(thisNode.getParent(), reads, assignedreads);
			}
		}catch(NullPointerException e){
			// TODO figure out why this is happening!
			//IO_utils.printLineErr("ERROR: NullPointer at node "+thisNode.getName()+" ("+reads.size()+" reads)");
		}
	}


	private HashSet<String> tidyUp(HashSet<String> reads, HashSet<String> assignedreads){
		// remove reads that have been assigned from the read2species map
		int readsAssignedHere = 0;
		//Iterator<String> it = _readAssignment.keySet().iterator();
		Iterator<String> it = assignedreads.iterator();
		while(it.hasNext()){
			String thisReadID = it.next();

			if(_reads2species.containsKey(thisReadID)){
				Iterator<String> tmpSpeciesToModify = _reads2species.get(thisReadID).iterator();
				while(tmpSpeciesToModify.hasNext()){
					String thisSpecies = tmpSpeciesToModify.next();

					// remove the read allocated to this species
					if(_species2reads.get(thisSpecies).contains(thisReadID))
						_species2reads.get(thisSpecies).remove(thisReadID);

					// reduce read count for this species by 1
					if(_species2readCounts.containsKey(thisSpecies))
						_species2readCounts.put(thisSpecies, _species2readCounts.get(thisSpecies)-1);
				}

				_reads2species.remove(thisReadID);
				readsAssignedHere ++;
			}

			//
			if(reads.contains(thisReadID))
				reads.remove(thisReadID);
		}
		if(_verbose == 2)
			System.err.println("\tDone, assigned "+readsAssignedHere+" reads here, "+_readAssignment.size()+" in total");
		return reads;
	}




	/**
	 * 
	 * @param node
	 */
	private void processAlignments(NCBITaxonomy_Node node){
		//System.out.println(node.getName()+": nChildren="+node.getChildren().size());

		Iterator<NCBITaxonomy_Node> childIterator = node.getChildren().iterator();
		while(childIterator.hasNext()){
			NCBITaxonomy_Node thisChild = childIterator.next();

			//System.out.print(node.getName()+": child:"+thisChild.getName());

			HashSet<String> leafNodes = thisChild.getAllLeafSpeciesForNode();

			//System.out.println(": leafNodeCount="+leafNodes.size());

			// check to see if this branch has any alignments
			boolean branchHasReads = false;
			Iterator<String> leafNodeIterator = leafNodes.iterator();
			while(leafNodeIterator.hasNext()){
				if(_species2reads.containsKey(leafNodeIterator.next())){
					branchHasReads = true;
					break;
				}
			}

			//System.out.println(thisChild.getName()+" has reads: "+branchHasReads);

			// if there are alignments for this branch, continue...
			if(branchHasReads){
				// get remaining reads
				Iterator<String> reads = _reads2species.keySet().iterator();
				while(reads.hasNext()){
					String thisRead = reads.next();

					int alignedSpeciesUnderCurrentParent = 0;
					Iterator<String> species = _reads2species.get(thisRead).iterator();
					while(species.hasNext()){
						if(leafNodes.contains(species.next()))
							alignedSpeciesUnderCurrentParent ++;
					}

					double frac = (alignedSpeciesUnderCurrentParent+0.0)/(_reads2species.get(thisRead).size()+0.0);

					//System.out.println(thisRead+"\t"+thisChild.getName()+" ("+thisChild.getLevel()+"): "+frac);

					if(frac >= _minReaadFrac){ // this read still works at this level
						_readAssignment.put(thisRead, thisChild);
						//System.out.println(thisRead+"\t"+thisChild.getName()+" ("+thisChild.getLevel()+"): "+frac);
					}
				}
				processAlignments(thisChild);
			}
		}
	}


	HashSet<String> _readsToIgnore = new HashSet<String>();
	private HashMap<String, NCBITaxonomy_Node> _readAssignment = new HashMap<String, NCBITaxonomy_Node>();
	private HashMap<NCBITaxonomy_Node, Integer> _nodeReadCount = new HashMap<NCBITaxonomy_Node, Integer>();



	private HashMap<String, HashSet<String>> _reads2species = new HashMap<String, HashSet<String>>();
	private HashMap<String, HashSet<String>> _species2reads = new HashMap<String, HashSet<String>>();
	private HashMap<String, Integer> _species2readCounts = new HashMap<String, Integer>();












	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("dir").hasArg().withDescription("Path to the directory containing the NCBI taxonomy").create("taxonomyPath"));
		options.addOption(OptionBuilder.withArgName(".unique.txt").hasArg().withDescription("Path to the exogenous alignments output by exceRpt").create("alignments"));
		options.addOption(OptionBuilder.withArgName("[min=0, max=1]").hasArg().withDescription("[optional] minimum fraction of reads that must be contained within a subtree before increasing resolution to that subtree [default: 0.95]").create("frac"));
		options.addOption(OptionBuilder.withArgName("integer").hasArg().withDescription("[optional] process alignments in batches of reads [default: "+DEFAULT_BATCH_SIZE+"]").create("batchSize"));
		options.addOption(OptionBuilder.withArgName("decimal").hasArg().withDescription("[optional] Rapid-run mode: do not resolve alignments for leaf-nodes with fewer than this % of total reads [default: 0.0] [suggested 0.001]").create("min"));
		options.addOption(OptionBuilder.withDescription("[optional] Perform tree search \'top down\'.  This is guaranteed to find the correct assignment for each read, but can be very slow for a large number of alignments").create("topDown"));
		options.addOption(OptionBuilder.withDescription("If the inout alignments are agains exogenous rRNA sequences in the RDP").create("rdp"));
		options.addOption(OptionBuilder.withArgName("integer").hasArg().withDescription("[optional] output taxa information for branches with at least this number of alignments [default: "+DEFAULT_MINREADS+"]").create("minReads"));
		options.addOption(OptionBuilder.withDescription("Print verbose status messages").create("v"));
		options.addOption(OptionBuilder.withDescription("Print REALLY verbose status messages").create("vv"));
		return options;
	}

	public static final int DEFAULT_MINREADS = 3;
	public static final int DEFAULT_BATCH_SIZE = 200000;
	public int _verbose = 0;

	public static void main(String[] args) throws IOException, ParseException {
		String alignmentsPath;
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/TCGA-EW-A1P8.normal.unmapped_nonIBC_TCGA-EW-A1P8_ExogenousGenomicAlignments.sorted.unique.txt";
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/Sample_Yale-IBC-298-normal.unmapped_ExogenousGenomicAlignments.sorted.unique.txt";
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/PG0004515-BLD.hs37d5.unmapped_ExogenousGenomicAlignments.sorted.unique.txt";
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/A806WMABXX.unmapped_ExogenousGenomicAlignments.sorted.unique.txt";
		alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/test.txt";

		//alignmentsPath = "/Users/robk/Downloads/TEST.txt";
		//args = new String[]{"--taxonomyPath","/Users/robk/WORK/YALE_offline/ANNOTATIONS/taxdump", "--alignments",alignmentsPath, "--batchSize","200000", "-v", "-min","0.001"};
		//args = new String[]{"--taxonomyPath","/Users/robk/WORK/YALE_offline/ANNOTATIONS/taxdump", "--alignments",alignmentsPath, "--batchSize","200000", "-v", "-topDown"};



		alignmentsPath = "/Users/robk/Downloads/ExogenousRibosomalAlignments.head";
		//args = new String[]{"--taxonomyPath","/Users/robk/WORK/YALE_offline/ANNOTATIONS/taxdump", "--alignments",alignmentsPath, "--rdp", "--minReads","1", "--batchSize","200000", "--v"};

		CommandLine cmdArgs = ExceRpt_Tools.parseArgs(args, getCmdLineOptions());
		if(cmdArgs.hasOption("taxonomyPath") && cmdArgs.hasOption("alignments")){

			int verbose = 0;
			if(cmdArgs.hasOption("vv")){
				IO_utils.printLineErr("Printing REALLY verbose status messages");
				verbose = 2;
			}else if(cmdArgs.hasOption("v")){
				IO_utils.printLineErr("Printing verbose status messages");
				verbose = 1;
			}

			boolean topDown = false;
			if(cmdArgs.hasOption("topDown"))
				topDown = true;

			int dataType = DATA_TYPE_GENOME;
			if(cmdArgs.hasOption("rdp")){
				dataType = DATA_TYPE_RRNA;
				//topDown = true;
			}

			int minReads = DEFAULT_MINREADS;
			if(cmdArgs.hasOption("minReads")){
				minReads = Integer.valueOf(cmdArgs.getOptionValue("minReads")).intValue();
			}

			double minPercentage = 0.0;
			if(cmdArgs.hasOption("min"))
				minPercentage = Double.valueOf(cmdArgs.getOptionValue("min")).doubleValue();

			ProcessExogenousAlignments engine = new ProcessExogenousAlignments(dataType, topDown, minPercentage/100.0, minReads, verbose);

			int batchSize = DEFAULT_BATCH_SIZE;
			if(cmdArgs.hasOption("batchSize"))
				batchSize = Integer.valueOf(cmdArgs.getOptionValue("batchSize")).intValue();

			double frac = 0.95;
			if(cmdArgs.hasOption("frac")){
				frac = Double.valueOf(cmdArgs.getOptionValue("frac")).doubleValue();
				engine.setMinReadFractionForTreeDescent(frac);
			}

			IO_utils.printLineErr("Summarising alignments in: \'"+cmdArgs.getOptionValue("alignments")+"\'");
			IO_utils.printLineErr("            in batches of: "+batchSize+" reads");
			IO_utils.printLineErr("        using taxonomy in: \'"+cmdArgs.getOptionValue("taxonomyPath")+"\', minFraction="+frac+", top-down:"+topDown+", for alignments above "+minPercentage+"% of total reads");

			engine.readTaxonomy(cmdArgs.getOptionValue("taxonomyPath"));

			engine.processAlignments(cmdArgs.getOptionValue("alignments"), batchSize);
		}
		else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(ExceRpt_Tools.ExceRpt_Tools_EXE_COMMAND+" ProcessExogenousAlignments", getCmdLineOptions());
			System.err.println();
		}
	}
}

package main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
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


	/**
	 * Main
	 * 
	 * @param taxonomyPath
	 * @param alignmentsPath
	 * @throws IOException
	 */
	public ProcessExogenousAlignments() throws IOException{  }
	
	
	/**
	 * 
	 * @param taxonomyPath
	 * @throws IOException
	 */
	public void readTaxonomy(String taxonomyPath) throws IOException{
		_taxonomy = new NCBITaxonomy_Engine(taxonomyPath);
	}
	
	
	/**
	 * 
	 * @param alignmentsPath
	 * @throws IOException
	 */
	public void processAlignments(String alignmentsPath) throws IOException{
		//IO_utils.printLineErr("Reading alignments");
		readAlignments(alignmentsPath);

		IO_utils.printLineErr("Processing alignments");
		processAlignments(_taxonomy.getRootNode());

		IO_utils.printLineErr("Assigned "+_readAssignment.size()+" of "+_reads2species.size()+" reads ("+Math.round((_readAssignment.size()+0.0)*1000.0/(_reads2species.size()+0.0))/10.0+"%) to a subset of the taxonomy");

		IO_utils.printLineErr("Quantifying");
		countReads();
		//IO_utils.printLineErr("Quantified taxonomy:");
		String indent = ">";
		int minReads = 3;
		_taxonomy.getRootNode().printSummary(indent, minReads);

		IO_utils.printLineErr("Done");
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



	private HashMap<String, ArrayList<String>> _reads2species = new HashMap<String, ArrayList<String>>();
	private HashMap<String, ArrayList<String>> _species2reads = new HashMap<String, ArrayList<String>>();
	private HashMap<String, Integer> _species2readCounts = new HashMap<String, Integer>();

	private void readAlignments(String alignmentsPath) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(alignmentsPath));
		String line;
		int ignoredReadCount = 0;
		HashSet<String> ignoredReads = new HashSet<String>();
		//HashSet<String> ignoredSpecies = new HashSet<String>();

		//int count = 0;
		while((line=br.readLine())!=null){
			//if(count < 10){
			//	System.out.println(line);
			//}

			String[] bits = line.split("\t");

			String readID = bits[0].trim();
			String species = bits[2].replace("_", " ").toLowerCase().trim();

			// check that this species is in the taxonomy
			//if(_taxonomy.getNodeName2nodeIndex().containsKey(species)){
			if(_taxonomy.containsNode(species)){

				if(!_reads2species.containsKey(readID))
					_reads2species.put(readID, new ArrayList<String>());
				_reads2species.get(readID).add(species);

				if(!_species2reads.containsKey(species)){
					_species2reads.put(species, new ArrayList<String>());
					_species2readCounts.put(species, 0);
				}
				_species2reads.get(species).add(readID);

				int thisCount = _species2readCounts.get(species)+1;
				_species2readCounts.put(species, thisCount);



				//System.out.println(species+": "+_species2readCounts.get(species));
			}else{
				//IO_utils.printLineErr("\'"+species+"\' not found in the taxonomy, ignoring");
				ignoredReads.add(readID);
				//ignoredSpecies.add(species);
			}

		}

		// count reads that are removed completely due to taxonomy ambiguity
		Iterator<String> it = ignoredReads.iterator();
		while(it.hasNext()){
			if(!_reads2species.containsKey(it.next()))
				ignoredReadCount ++;
		}

		IO_utils.printLineErr("Ignored "+ignoredReads.size()+" alignments to species that could not be found in the taxonomy ("+ignoredReadCount+" reads removed entirely)");
		IO_utils.printLineErr("Detected "+_species2reads.size()+" possible species from "+_reads2species.size()+" aligned reads");
		//System.out.println(_maxCountSpecies+": "+_maxCount);
		br.close();
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("dir").hasArg().withDescription("Path to the directory containing the NCBI taxonomy").create("taxonomyPath"));
		options.addOption(OptionBuilder.withArgName(".unique.txt").hasArg().withDescription("Path to the exogenous alignments output by exceRpt").create("alignments"));
		options.addOption(OptionBuilder.withArgName("[min=0, max=1]").hasArg().withDescription("[optional] minimum fraction of reads that must be contained within a subtree before increasing resolution to that subtree").create("frac"));
		return options;
	}


	public static void main(String[] args) throws IOException, ParseException {
		//String alignmentsPath;
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/TCGA-EW-A1P8.normal.unmapped_nonIBC_TCGA-EW-A1P8_ExogenousGenomicAlignments.sorted.unique.txt";
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/Sample_Yale-IBC-298-normal.unmapped_ExogenousGenomicAlignments.sorted.unique.txt";
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/PG0004515-BLD.hs37d5.unmapped_ExogenousGenomicAlignments.sorted.unique.txt";
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/A806WMABXX.unmapped_ExogenousGenomicAlignments.sorted.unique.txt";
		//alignmentsPath = "/Users/robk/WORK/YALE_offline/tmp/Lajos/exogenousTEST/test.txt";

		//args = new String[]{"--taxonomyPath","/Users/robk/WORK/YALE_offline/ANNOTATIONS/taxdump", "--alignments",alignmentsPath};
		//String taxonomyPath = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/taxdump";
		

		CommandLine cmdArgs = ExceRpt_Tools.parseArgs(args, getCmdLineOptions());
		if(cmdArgs.hasOption("taxonomyPath") && cmdArgs.hasOption("alignments")){
			
			ProcessExogenousAlignments engine = new ProcessExogenousAlignments();

			double frac = 0.95;
			if(cmdArgs.hasOption("frac")){
				frac = Double.valueOf(cmdArgs.getOptionValue("frac")).doubleValue();
				engine.setMinReadFractionForTreeDescent(frac);
			}
			
			IO_utils.printLineErr("Summarising alignments in \'"+cmdArgs.getOptionValue("alignments")+"\'");
			IO_utils.printLineErr("        using taxonomy in \'"+cmdArgs.getOptionValue("taxonomyPath")+"\' and minFraction="+frac);
			
			engine.readTaxonomy(cmdArgs.getOptionValue("taxonomyPath"));
			
			engine.processAlignments(cmdArgs.getOptionValue("alignments"));
		}
		else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(ExceRpt_Tools.ExceRpt_Tools_EXE_COMMAND+" ProcessExogenousAlignments", getCmdLineOptions());
			System.err.println();
		}
	}
}

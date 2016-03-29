import java.util.ArrayList;

/**
 * 
 * @author Nick Fisk
 * A node containing information about the 
 * organism used in searching a graph to recover a tree.
 * currently used only for a UPGMA approach, but
 * is extendable to other approaches. 
 */
public class Node {
	String newick;
	ArrayList<String> names=new ArrayList<String>();
	
	public Node(String name){
		this.names.add(name);
		name=name.replaceAll("[()]", "");
		name=name.replaceAll(">", "");
		//name=name.replaceAll("\\","");
		//name=name.replaceAll("/","");
		if(name.length()>30){
			name=name.substring(0, 29);
		}
		name=name.replaceAll(" ", "_");
		name=name.replaceAll("-", "_");
		name.replaceAll(",", "_");
		this.newick=name;
	}
	

	
	/**
	 * Constructs a node from two other nodes and a distance
	 * @param node1, a node to be joined to node2
	 * @param node2, a node to be joined to node1
	 * @param dist the distance betwwen the two input nodes
	 * 
	 */
	@SuppressWarnings("unchecked")
	public Node(Node node1, Node node2, double dist){
		this.newick="("+node1.newick+":"+dist/2+","+node2.newick+":"+dist/2+")";
		this.names=new ArrayList<String>(node1.names);
		this.names.addAll((ArrayList<String>)node2.names.clone());
	}
	/**
	 * Find the distance between any two nodes. For nodes that 
	 * consist of many nodes, the overall average of all the nodes is used,
	 * not just the average of the two nodes being compared
	 * @param n1, the first node or cluster of nodes
	 * @param n2, the second node or cluster of nodes
	 * @return the distance between the two nodes
	 */
	public static double calcDistance(Node n1, Node n2){
		ArrayList<String>names1=new ArrayList<String>(n1.names);
		ArrayList<String>names2=new ArrayList<String>(n2.names);
		double score=0.0;
		int count=0;
		for(String i: names1){
			for(String j: names2){
				score+=ScrawkovPHY.initScores.get(i).get(j);
				count+=1;
			}
		}
		return score/count;
	}
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString(){
		return(String.valueOf(this.names));
	}
	
	@Override
	public boolean equals(Object o){
		if(o==this){
			return true;
		}
		if(!(o instanceof Node)){
			return false;
		}
		Node o2=(Node)o;
		for(String n: this.names){
			if(o2.names.contains(n)==false){
				return false;
			}
		}
		for(String n: o2.names){
			if(this.names.contains(n)==false){
				return false;
			}
		}
		return true;
	}
	/**
	 * needed for the equals override
	 */
	public int hashCode(){
		return this.names.hashCode();
	}
}

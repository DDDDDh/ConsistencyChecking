package cn.edu.nju.moon.consistency.ui;

import java.io.File;

import cn.edu.nju.moon.consistency.checker.ReadIncChecker;
import cn.edu.nju.moon.consistency.model.observation.RawObservation;
import cn.edu.nju.moon.consistency.model.operation.BasicOperation;
import cn.edu.nju.moon.consistency.model.process.RawProcess;

/**
 * @description UI for {@link ReadIncChecker}; visualization of the data structures 
 * 	used in checking algorithm
 * 
 *  The UI is based on GraphViz (with .dot).
 *  
 * @author hengxin
 * @date 2013-1-2
 * @design_pattern Singleton
 */
public class DotUI
{
	private GraphViz viz = null;
	private final String type = "pdf";
	private File out = null;
	
	private static DotUI instance = null;
	private DotUI()
	{
		this.viz = new GraphViz();
		this.viz.addln(viz.start_graph());
	}
	
	public static DotUI getInstance()
	{
		if (instance == null)
			instance = new DotUI();
		
		return instance;
	}
	
	/**
	 * execute external GraphViz program
	 * @param out file for output
	 */
	public void execute(String out)
	{
		this.viz.addln(this.viz.end_graph());
    	
    	// for test
    	System.out.println(viz.getDotSource());
		
		this.out = new File("data/" + out + "." + this.type);
		System.out.println(this.out.getAbsolutePath());
		viz.writeGraphToFile(viz.getGraph(viz.getDotSource(), type), this.out);
	}
	
	/**
	 * visualize {@link RawObservation}
	 * 
	 * @param riob {@link RawObservation} to visualize
	 */
    public void visual_ob(RawObservation ob)
    {
    	int size = ob.getSize();
    	
    	viz.addln("ranksep = 1.0; size = \"10,10\";");	// drawing with constrained ranks
    	viz.addln("{");	// pids
    		viz.addln("node [shape = plaintext, fontsize = 20];");
    		for (int i = 0; i < size - 1; i++)
    			viz.add(i + " -> ");
    		viz.addln((size - 1) + ";");
    	viz.addln("}");
    	
    	viz.addln("node [shape = box];");	// shape of nodes
    	
    	for (int pid : ob.getProcMap().keySet())
    	{
    		this.visual_proc(ob.getProcess(pid));
    	}
    }
    
    /**
     * visualize {@link RawProcess}
     * 
     * @param ripro {@link RawProcess} to visualize
     */
    private void visual_proc(RawProcess ripro)
    {
    	viz.add("{");
    		viz.add("rank = same;");	// in the same process
    		viz.add(ripro.getPid() + ";");
    		for (BasicOperation op : ripro.getOpList())
    			this.visual_op(op);
    	viz.addln("}");
    }
    
    /**
     * visualize {@link BasicOperation}
     * @param op {@link BasicOperation} to visualize
     */
    private void visual_op(BasicOperation op)
    {
    	viz.add(op.toString() + ";");
    }
    
    /**
     * add program order edge from @param from_op to @param to_op
     * 
     * @param from_op source of edge
     * @param to_op   target of edge
     */
    public void addPOEdge(BasicOperation from_op, BasicOperation to_op)
    {
//    	this.viz.addComment("Program Order:");
    	this.viz.addln(from_op.toString() + " -> " + to_op.toString() + ";");
    }
    
    /**
     * add WriteTo order edge from @param from_op to @param to_op
     * @param from_op source of edge
     * @param to_op   target of edge
     */
    public void addWritetoEdge(BasicOperation from_op, BasicOperation to_op)
    {
//    	this.viz.addComment("Write to Order:");
    	this.viz.addln(from_op.toString() + " -> " + to_op.toString() + "[color = blue];");
    }
    
    public void addWprimeWREdge(BasicOperation from_op, BasicOperation to_op)
    {
    	this.viz.addComment("W'WR Order:");
    	this.viz.addln(from_op.toString() + " -> " + to_op.toString() + "[style = dashed, color = red];");
    }
}
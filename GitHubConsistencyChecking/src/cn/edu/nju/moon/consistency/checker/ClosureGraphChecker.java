package src.cn.edu.nju.moon.consistency.checker;

import src.cn.edu.nju.moon.consistency.model.observation.ClosureObservation;
import src.cn.edu.nju.moon.consistency.model.observation.BasicObservation;
import src.cn.edu.nju.moon.consistency.schedule.ISchedule;
import src.cn.edu.nju.moon.consistency.ui.DotUI;

/**
 * @description OperationGraphChecker is a PRAM Consistency checking algorithm
 *   which based on ``closure'' idea
 *   
 * @author hengxin
 * @date 2013-1-8
 */
public class ClosureGraphChecker extends Checker
{

	public ClosureGraphChecker(BasicObservation rob)
	{
		super(rob);
	}

	/**
	 * Constructor
	 * @param rob	{@link BasicObservation} to check
	 * @param name	name for visualization in {@link DotUI}
	 */
	public ClosureGraphChecker(BasicObservation rob, String name)
	{
		super(rob, name);
	}
	
	/**
	 * Constructor
	 * @param bob	{@link BasicObservation} to check
	 * @param name	name for visualization in {@link DotUI}
	 * @param s {@link ISchedule}: record for the checking result
	 */
	public ClosureGraphChecker(BasicObservation bob, String name, ISchedule s)
	{
		super(bob, name, s);
	}
	
	/**
	 * get the appropriate {@link ClosureObservation} with respect to @param masterPid
	 */
	@Override
	protected BasicObservation getMasterObservation(int masterPid)
	{
//		System.out.println("oh?");
		return new ClosureObservation(masterPid, super.rob);
	}

	/**
	 * @description operation-graph checking algorithm (based on iterative closure)
	 * 
	 * @param rob observation to check (i.e., {@link ClosureObservation})
	 */
	@Override
	protected boolean check_part(BasicObservation rob)
	{
		boolean consistent = true;
		
		ClosureObservation clob = (ClosureObservation) rob;

		/**
		 * @modified hengxin on 2013-1-9
		 * @modification remove "clob.preprocessing()" here
		 * @reason it has been called by super class {@link Checker}#check() (as a Template Method)
		 */
//		// 1) add program order edge and read-write mapping
//		clob.preprocessing();

		/**
		 *  2) iteration of applying two rules until no new edges are added,
		 *  in other words, a fixed point is reached
		 */
		this.loopTime = 0;

		boolean changed_wprime2w = false;
		do
		{
			this.loopTime++;
			clob.closure();
			changed_wprime2w = clob.apply_wprime2w_edge_rule();
			if (clob.cycle_check())	// cycle detection
			{
				consistent = false;
				break;
			}
		}while(changed_wprime2w);
		System.out.println("Loop time of closure alg for process" + clob.getMasterPid() + ": " + this.loopTime);

//		if (consistent)	// satisfying PRAM Consistency; construct a witness
//			schedule.constructExecution(master_cls_ob);

		// ui
		DotUI.getInstance().execute("operationgraph/" + name + "_" + clob.getMasterPid());
		
		return consistent;
	}

}

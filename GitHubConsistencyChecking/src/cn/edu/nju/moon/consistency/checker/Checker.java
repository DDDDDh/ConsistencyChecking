package src.cn.edu.nju.moon.consistency.checker;

import org.apache.commons.lang3.RandomStringUtils;

import src.cn.edu.nju.moon.consistency.model.observation.BasicObservation;
import src.cn.edu.nju.moon.consistency.model.observation.ReadIncObservation;
import src.cn.edu.nju.moon.consistency.model.process.BasicProcess;
import src.cn.edu.nju.moon.consistency.schedule.ISchedule;
import src.cn.edu.nju.moon.consistency.schedule.WeakSchedule;
import src.cn.edu.nju.moon.consistency.ui.DotUI;

import src.cn.edu.nju.moon.consistency.model.operation.*;

/**
 * @description Consistency checking algorithm is responsible for implementing
 * 	{@link IChecker}
 * 
 * @author hengxin
 * @date 2012-12-8
 * 
 * @modified hengxin on 2013-1-5
 * @reason using Template Method design pattern
 */
public abstract class Checker
{
	protected BasicObservation rob = null;	/** {@link BasicObservation} to check **/
	protected String name = "";		/** for {@link DotUI}; the name of file for visualization **/
	private ISchedule schedule = null;	/** record of the checking result */
	public int loopTime = 0;
	public int skipLoop = 0;
	/**
	 * Constructor
	 * @param bob {@link BasicObservation} to check
	 */
	public Checker(BasicObservation bob)
	{
		this.rob = bob;
		this.name = RandomStringUtils.random(8);
	}
	
	/**
	 * Constructor
	 * @param riob	{@link BasicObservation} to check
	 * @param name	for {@link DotUI}; the name of file for visualization
	 */
	public Checker(BasicObservation bob, String name)
	{
		this.rob = bob;
		this.name = name;
	}
	
	/**
	 * Constructor
	 * @param bob {@link BasicObservation} to check
	 * @param s {@link ISchedule}: record for the checking result
	 */
	public Checker(BasicObservation bob, ISchedule s)
	{
		this(bob);
		this.schedule = s;
	}
	
	/**
	 * Constructor
	 * @param bob {@link BasicObservation} to check
	 * @param name for {@link DotUI}; the name of file for visualization
	 * @param s  {@link ISchedule}: record for the checking result
	 */
	public Checker(BasicObservation bob, String name, ISchedule s)
	{
		this(bob, name);
		this.schedule = s;
	}
	
	/**
	 * check whether {@link #rob} satisfies PRAM Consistency
	 * @return true, if {@link #rob} satisfies PRAM Consistency; false, otherwise.
	 * 
	 * Template Method design pattern
	 */
	public final boolean check()
	{
		boolean consistent = true;
		
		int pids = this.rob.getProcNum();
		BasicObservation mob = null;

		for (int pid = 0; pid < pids; pid++)
		{
			boolean partial_consistent = true;
			System.out.println("dealing with process" + pid);
			mob = this.getMasterObservation(pid);
			if (this.trivial_check(mob))	/** pass the simple check */
			{
//				System.out.println("Pass trivial check!");

				if (! this.check_part(mob))	/** process with pid does not satisfy consistency condition **/
				{
					System.out.println("detect cycle at process " + pid);
					consistent = false;
					partial_consistent = false;
				}
				System.out.println("Loop time for process" + pid +":" + this.loopTime);
				System.out.println("Skip time for process" + pid +":" + this.skipLoop);
				this.loopTime = 0;
				this.skipLoop = 0;

			}
			else
			{
				consistent = false;
				partial_consistent = false;
			}
			
			if (this.schedule != null && partial_consistent)
				this.schedule.constructView(mob);	/** the process with pid satisfies consistency condition; construct a view for it */
		}
		
		return consistent;
	}
	
	/**
	 * @return {@link #schedule}
	 */
	public ISchedule getSchedule()
	{
		return this.schedule;
	}
	
	/**
	 * simple check for PRAM Consistency
	 * 1. nullCheck: no operations in the process to be checked; 
	 * 		it is trivially PRAM Consistent
	 * 2. readLaterWrite: some READ reads later WRITE in the same process; 
	 * 		it does not satisfy PRAM Consistency
	 * 
	 * @param ob {@link BasicObservation} to check
	 * @return true, if it passes the simple check; false, otherwise.
	 */
	private boolean trivial_check(BasicObservation ob)
	{
		if (ob.nullCheck())	/** no operations in the process to be checked; it is trivially PRAM Consistent **/
		{
//			System.out.println("Null Check: true");
			return true;
		}

		
		ob.preprocessing();	// preprocessing: program order and write to order

		if (ob.readLaterWrite())	/** some READ reads later WRITE in the same process; it does not satisfy PRAM Consistency **/
		{
//			System.err.println("Read late write: false");
			return false;
		}
		
		return true;
	}
	
	/**
	 * @param masterPid pid of {@link BasicProcess} to check against PRAM Consistency
	 * @return specific subclass of {@link BasicObservation} with @param masterPid to check 
	 */
	protected abstract BasicObservation getMasterObservation(int masterPid);
	
	/**
	 * check with respect to some process against PRAM Consistency
	 * 
	 * @return true, if this process satisfies PRAM Consistency; false, otherwise. 
	 */
	protected abstract boolean check_part(BasicObservation masterObservation);
}

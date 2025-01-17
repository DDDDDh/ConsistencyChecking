package src.cn.edu.nju.moon.consistency.model.process;

import java.util.ArrayList;
import java.util.List;

import src.cn.edu.nju.moon.consistency.model.observation.BasicObservation;
import src.cn.edu.nju.moon.consistency.model.operation.BasicOperation;
import src.cn.edu.nju.moon.consistency.model.operation.ClosureOperation;
import src.cn.edu.nju.moon.consistency.model.operation.RawOperation;
import src.cn.edu.nju.moon.consistency.model.operation.ReadIncOperation;
import src.cn.edu.nju.moon.consistency.model.operation.factory.IOperationTransformer;

/**
 * @author hengxin
 * @date 2012-12-6
 * 
 * @description raw process is a list of {@link BasicOperation}
 *
 * @modified hengxin on 2013-1-8
 * @reason refactor: no RawProcess (in contrast to {@link RawOperation} 
 * 	change to {@link BasicProcess}
 */
public class BasicProcess
{
	protected int pid = -1;
	protected List<BasicOperation> opList = new ArrayList<BasicOperation>();
	
	public BasicProcess(int pid)
	{
		this.pid = pid;
	}
	
	/**
	 * add the {@link BasicOperation} to the {@link BasicProcess}
	 * and set its pid to the id of {@link BasicProcess}
	 * 
	 * @param op {@link BasicOperation} to be added
	 */
	public void addOperation(BasicOperation op)
	{
		op.setPid(pid);
		this.opList.add(op);
	}
	
	/**
	 * @return {@link #pid}
	 */
	public int getPid()
	{
		return this.pid;
	}
	
	/**
	 * return the {@link BasicOperation} with index
	 * @param index index of {@link BasicOperation} to return
	 */
	public BasicOperation getOperation(int index)
	{
		return this.opList.get(index);
	}
	
	/**
	 * @return the immutable field {@link #opList}
	 * @see {@link #getOpListCopy()} 
	 */
	public final List<BasicOperation> getOpList()
	{
		return this.opList;
	}
	
	/**
	 * @return the mutable copy of field {@link #opList}
	 * @see {@link #getOpList()}
	 */
	public List<BasicOperation> getOpListCopy()
	{
		return new ArrayList<BasicOperation>(this.opList);
	}
	
	/**
	 * @return number of {@link BasicOperation}s in this {@link BasicProcess}
	 */
	public int getOpNum()
	{
		return this.opList.size();
	}
	
	/**
	 * establish "program order" between {@link ClosureOperation}s 
	 * in the same {@link ClosureProcess}
	 */
	public void establishProgramOrder()
	{
		if (this.opList.size() == 0)
			return;
		
		BasicOperation preOp = this.opList.get(0);
		BasicOperation curOp = null;
		int size = this.opList.size();
		
		for (int index = 1; index < size; index++)
		{
			curOp = this.opList.get(index);
			preOp.setProgramOrder(curOp);
			preOp = curOp;
		}
	}
	
	/**
	 * establish WritetoOrder: D(R) => R
	 * @param bob {@link BasicObservation} from which you can retrieve any WRITE
	 */
	public void establishWritetoOrder(BasicObservation bob)
	{
		List<BasicOperation> opList = this.opList;
		BasicOperation rclop = null;
		BasicOperation wclop = null;
		int size = opList.size();
		for (int index = 0; index < size; index++)
		{	
			rclop = opList.get(index);
			if(rclop.isReadOp())	
			{
				wclop = bob.getDictatingWrite(rclop);
				wclop.addWritetoOrder(rclop);
			}
		}
	}
	
	/**
	 * does some READ {@link BasicOperation} read value from later WRITE
	 * {@link BasicOperation} on the same {@link BasicProcess}
	 * 
	 * @return true, if it does; false, otherwise.
	 */
	public boolean readLaterWrite()
	{
		BasicOperation wriop = null;
		for (BasicOperation riop : this.opList)
		{
			if (riop.isReadOp())	// check every READ
			{
				wriop = riop.getReadfromWrite();
				if (riop.getPid() == wriop.getPid() && riop.getIndex() < wriop.getIndex())
					return true;
			}
		}
		return false;
	}
	
	/**
	 * filter <code>this</code> {@link BasicProcess}, 
	 * transform the remaining operations, and
	 * fill them into other process: @param proc
	 *   
	 * @param masterPid process with masterPid will keep its READ operations
	 * @param proc process to be constructed (filled)
	 * @param op_trans used to transform operations into appropriate type
	 */
	public void filter_fill(int masterPid, BasicProcess proc, IOperationTransformer op_trans)
	{
		for (BasicOperation bop : this.getOpListCopy())
		{	
			BasicOperation op = null; 
			
			if (bop.isWriteOp() || (bop.isReadOp() && proc.getPid() == masterPid) )
			{
				op = op_trans.transform(bop);
				op.setIndex(proc.getOpNum());
				proc.addOperation(op);
			}
		}
	}

	public void filter_fill_causal(int masterPid, BasicProcess proc, IOperationTransformer op_trans){ //与上面函数的区别在于保留了非主线程上的读操作
		for(BasicOperation bop: this.getOpListCopy()){
			BasicOperation op = null;
			op = op_trans.transform(bop);
			op.setIndex(proc.getOpNum());
			proc.addOperation(op);
		}
	}
	
}

# GraphOps
A Dataflow library for graph analytics acceleration

## Enumeration of GraphOps Components

### Data-Handling Blocks
Data blocks are the primary GraphOps components.  They handle incoming data
streams, perform arithmetic operations, and route outputs to memory or
subsequent blocks:

* **ForAllPropRdr** issues memory requests for all neighbor property sets
in the graph.  In order to do this, it first reads all the row pointers in the graph.
The incoming row pointer data are used to issue individual memory requests for
each set of neighbor properties.  Metadata about the requested neighbor
properties are emitted as an output to be processed by the subsequent block.
* **NbrPropRed** performs a reduction on a vertex's neighbor set.
The unit receives the neighbor property data as a data stream from memory.
Each set of neighbor properties is accompanied by a metadata packet as
an input to the kernel from a preceding block (e.g. ForAllPropRdr).  The
metadata is used to consume the correct amount of data from the incoming
neighbor property data stream.  For each neighbor property set, an accumulating
reduction is performed on the data and the result is emitted as an output along
with accompanying metadata.
* **ElemUpdate** is used to update property values in the graph data
structure.  The unit receives a vertex reference and an updated value as input.
It issues memory read requests for the requisite memory locations and memory
update requests for the updated values.
* **AllNodePropRed** reads property values for the entire graph and
performs actions based upon whether the values satisfy a condition.  Vertices
whose properties satisfy the condition are emitted to be processed in subsequent
blocks.  The properties themselves may also be optionally used for
computation (e.g. reduction) within this block.
* **NbrPropRdr** requests the properties of the neighbor set for
a given vertex, which is supplied as input.  This is the single-vertex version
of ForAllPropRdr.
* **SetReader** reads a set of vertices from memory and emits them as an
output stream to subsequent blocks.  This block is useful for algorithms
that generate intermediate working sets, e.g. *frontier sets*, between
iterations.
* **SetWriter** accumulates a working set of vertices and streams them out
to memory.  As with SetReader, this functionality is useful when creating
working sets.
* **NbrPropFilter** issues memory requests for the properties of a neighbor
set.  It filters the property values according to some condition and emits
properties which satisfy the condition, along with accompanying metadata.
* **GlobNbrRed** is used to perform a global reduction across an entire
subset of the graph.  The unit takes in a stream of property values, along with
accompanying metadata.  It uses the metadata to filter out non-applicable
properties.  For the duration of the execution/iteration, it performs an
accumulating reduction on the data.
* **VertexReader** issues memory requests for row pointers for a single
vertex and emits metadata.
* **NbrSetReader** uses row pointers and metadata to issue memory requests
for one vertex's list of neighbors.

### Control Blocks
In GraphOps, the majority of the logic is amenable to dataflow.  One key
reason for this is that feedback control is rare.  There are situations,
however, that call for more intricate control difficult to express
without state machines.  The control blocks are embedded inside the data
blocks and are responsible for handling these situations.  They are as
follows:
* **QRdrPktCntSM** handles control logic for input buffers in the
data blocks.  A common use case occurs in the following situation:  A
metadata input datum dictates how many memory packets belong to a given neighbor
set.  In this case, the QRdrPktCntSM block handles the counting of packets on
the memory data input and instructs the data block when to move on to the
next neighbor set.  This unit is also used for flow control, emitting a
*stall* output signal when metadata input buffers are nearing capacity.
* **UpdQueueSM** handles control logic for updating a graph property for
all nodes.  This unit assumes that the properties are being updated sequentially
and makes use of heavy coalescing to minimize the number of update requests sent
to memory.
* **CoalesceSM** also handles logic for updating a graph property.
However, this version does not assume in-order vertex updates and thus does not
coalesce as efficiently.  Best-effort coalescing buffer logic is built into the
control unit.
* **FifoKernelSM** is a control wrapper for a standard FIFO block.
It provides an additional *dataReady* control signal that is necessary to
construct more sophisticated queuing structures in the data blocks.
* **MemUnitSM** handles requests involving very large data sizes.  The
hardware platform underlying the GraphOps system may have a maximum limit for
size of request, so control logic is needed to issue multiple requests in this
case.  The unit includes input buffering to prevent subsequent requests from
being dropped while a large request is being issued.


### Utility Blocks
Additional logic is needed to properly interface with the memory system
and the host platform.  These are realized via the utility blocks:
* **EndSignal** monitors *done* signals for all data blocks and
issues a special interrupt request to halt execution when all units are
finished.
* **MemUnit** provides a simplified memory interface to the data
blocks.  It compiles memory profiling information, watches for
end-of-execution interrupt requests, and includes control logic for handling
very large memory requests.


## License

GraphOps is licensed under the [MIT License](http://opensource.org/licenses/MIT).

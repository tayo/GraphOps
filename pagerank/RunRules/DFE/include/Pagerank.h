/**\file */
#ifndef SLIC_DECLARATIONS_Pagerank_H
#define SLIC_DECLARATIONS_Pagerank_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface readLMem ----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'readLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_tocpu The stream should be of size (param_size * 4) bytes.
 */
void Pagerank_readLMem(
	uint64_t param_size,
	uint64_t param_start,
	uint32_t *outstream_tocpu);

/**
 * \brief Basic static non-blocking function for the interface 'readLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_tocpu The stream should be of size (param_size * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Pagerank_readLMem_nonblock(
	uint64_t param_size,
	uint64_t param_start,
	uint32_t *outstream_tocpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'readLMem'
 * 
 */
typedef struct { 
	uint64_t param_size; /**<  [in] Interface Parameter "size". */
	uint64_t param_start; /**<  [in] Interface Parameter "start". */
	uint32_t *outstream_tocpu; /**<  [out] The stream should be of size (param_size * 4) bytes. */
} Pagerank_readLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Pagerank_readLMem_run(
	max_engine_t *engine,
	Pagerank_readLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'readLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_readLMem_run_nonblock(
	max_engine_t *engine,
	Pagerank_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Pagerank_readLMem_run_group(max_group_t *group, Pagerank_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_readLMem_run_group_nonblock(max_group_t *group, Pagerank_readLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Pagerank_readLMem_run_array(max_engarray_t *engarray, Pagerank_readLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_readLMem_run_array_nonblock(max_engarray_t *engarray, Pagerank_readLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Pagerank_readLMem_convert(max_file_t *maxfile, Pagerank_readLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*--------------------------- Interface writeLMem ----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'writeLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size (param_size * 4) bytes.
 */
void Pagerank_writeLMem(
	uint64_t param_size,
	uint64_t param_start,
	const uint32_t *instream_fromcpu);

/**
 * \brief Basic static non-blocking function for the interface 'writeLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size (param_size * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Pagerank_writeLMem_nonblock(
	uint64_t param_size,
	uint64_t param_start,
	const uint32_t *instream_fromcpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'writeLMem'
 * 
 */
typedef struct { 
	uint64_t param_size; /**<  [in] Interface Parameter "size". */
	uint64_t param_start; /**<  [in] Interface Parameter "start". */
	const uint32_t *instream_fromcpu; /**<  [in] The stream should be of size (param_size * 4) bytes. */
} Pagerank_writeLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Pagerank_writeLMem_run(
	max_engine_t *engine,
	Pagerank_writeLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'writeLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_writeLMem_run_nonblock(
	max_engine_t *engine,
	Pagerank_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Pagerank_writeLMem_run_group(max_group_t *group, Pagerank_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_writeLMem_run_group_nonblock(max_group_t *group, Pagerank_writeLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Pagerank_writeLMem_run_array(max_engarray_t *engarray, Pagerank_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_writeLMem_run_array_nonblock(max_engarray_t *engarray, Pagerank_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Pagerank_writeLMem_convert(max_file_t *maxfile, Pagerank_writeLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] param_NumNodes Interface Parameter "NumNodes".
 * \param [in] param_StartCnt Interface Parameter "StartCnt".
 * \param [in] param_StopCnt Interface Parameter "StopCnt".
 * \param [in] param_d Interface Parameter "d".
 * \param [in] param_nodeAddr Interface Parameter "nodeAddr".
 * \param [in] param_prTerm Interface Parameter "prTerm".
 * \param [in] param_propAddr Interface Parameter "propAddr".
 * \param [in] param_repAddr Interface Parameter "repAddr".
 * \param [in] param_uDVal Interface Parameter "uDVal".
 * \param [out] outscalar_MemUnit0_numBrsts Output scalar parameter "MemUnit0.numBrsts".
 * \param [out] outscalar_MemUnit1_numBrsts Output scalar parameter "MemUnit1.numBrsts".
 * \param [out] outscalar_MemUnit2_numBrsts Output scalar parameter "MemUnit2.numBrsts".
 * \param [out] outscalar_MemUnit3_numBrsts Output scalar parameter "MemUnit3.numBrsts".
 * \param [out] outscalar_MemUnit4_numBrsts Output scalar parameter "MemUnit4.numBrsts".
 */
void Pagerank(
	uint32_t param_NumNodes,
	uint64_t param_StartCnt,
	uint64_t param_StopCnt,
	float param_d,
	uint32_t param_nodeAddr,
	float param_prTerm,
	uint32_t param_propAddr,
	uint32_t param_repAddr,
	uint32_t param_uDVal,
	uint64_t *outscalar_MemUnit0_numBrsts,
	uint64_t *outscalar_MemUnit1_numBrsts,
	uint64_t *outscalar_MemUnit2_numBrsts,
	uint64_t *outscalar_MemUnit3_numBrsts,
	uint64_t *outscalar_MemUnit4_numBrsts);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_NumNodes Interface Parameter "NumNodes".
 * \param [in] param_StartCnt Interface Parameter "StartCnt".
 * \param [in] param_StopCnt Interface Parameter "StopCnt".
 * \param [in] param_d Interface Parameter "d".
 * \param [in] param_nodeAddr Interface Parameter "nodeAddr".
 * \param [in] param_prTerm Interface Parameter "prTerm".
 * \param [in] param_propAddr Interface Parameter "propAddr".
 * \param [in] param_repAddr Interface Parameter "repAddr".
 * \param [in] param_uDVal Interface Parameter "uDVal".
 * \param [out] outscalar_MemUnit0_numBrsts Output scalar parameter "MemUnit0.numBrsts".
 * \param [out] outscalar_MemUnit1_numBrsts Output scalar parameter "MemUnit1.numBrsts".
 * \param [out] outscalar_MemUnit2_numBrsts Output scalar parameter "MemUnit2.numBrsts".
 * \param [out] outscalar_MemUnit3_numBrsts Output scalar parameter "MemUnit3.numBrsts".
 * \param [out] outscalar_MemUnit4_numBrsts Output scalar parameter "MemUnit4.numBrsts".
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Pagerank_nonblock(
	uint32_t param_NumNodes,
	uint64_t param_StartCnt,
	uint64_t param_StopCnt,
	float param_d,
	uint32_t param_nodeAddr,
	float param_prTerm,
	uint32_t param_propAddr,
	uint32_t param_repAddr,
	uint32_t param_uDVal,
	uint64_t *outscalar_MemUnit0_numBrsts,
	uint64_t *outscalar_MemUnit1_numBrsts,
	uint64_t *outscalar_MemUnit2_numBrsts,
	uint64_t *outscalar_MemUnit3_numBrsts,
	uint64_t *outscalar_MemUnit4_numBrsts);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint32_t param_NumNodes; /**<  [in] Interface Parameter "NumNodes". */
	uint64_t param_StartCnt; /**<  [in] Interface Parameter "StartCnt". */
	uint64_t param_StopCnt; /**<  [in] Interface Parameter "StopCnt". */
	float param_d; /**<  [in] Interface Parameter "d". */
	uint32_t param_nodeAddr; /**<  [in] Interface Parameter "nodeAddr". */
	float param_prTerm; /**<  [in] Interface Parameter "prTerm". */
	uint32_t param_propAddr; /**<  [in] Interface Parameter "propAddr". */
	uint32_t param_repAddr; /**<  [in] Interface Parameter "repAddr". */
	uint32_t param_uDVal; /**<  [in] Interface Parameter "uDVal". */
	uint64_t *outscalar_MemUnit0_numBrsts; /**<  [out] Output scalar parameter "MemUnit0.numBrsts". */
	uint64_t *outscalar_MemUnit1_numBrsts; /**<  [out] Output scalar parameter "MemUnit1.numBrsts". */
	uint64_t *outscalar_MemUnit2_numBrsts; /**<  [out] Output scalar parameter "MemUnit2.numBrsts". */
	uint64_t *outscalar_MemUnit3_numBrsts; /**<  [out] Output scalar parameter "MemUnit3.numBrsts". */
	uint64_t *outscalar_MemUnit4_numBrsts; /**<  [out] Output scalar parameter "MemUnit4.numBrsts". */
} Pagerank_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Pagerank_run(
	max_engine_t *engine,
	Pagerank_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_run_nonblock(
	max_engine_t *engine,
	Pagerank_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Pagerank_run_group(max_group_t *group, Pagerank_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_run_group_nonblock(max_group_t *group, Pagerank_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Pagerank_run_array(max_engarray_t *engarray, Pagerank_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Pagerank_run_array_nonblock(max_engarray_t *engarray, Pagerank_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Pagerank_convert(max_file_t *maxfile, Pagerank_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* Pagerank_init(void);

/* Error handling functions */
int Pagerank_has_errors(void);
const char* Pagerank_get_errors(void);
void Pagerank_clear_errors(void);
/* Free statically allocated maxfile data */
void Pagerank_free(void);
/* These are dummy functions for hardware builds. */
int Pagerank_simulator_start(void);
int Pagerank_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_Pagerank_H */


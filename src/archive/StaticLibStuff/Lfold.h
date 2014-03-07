#ifndef __VIENNA_RNA_PACKAGE_LFOLD_H__
#define __VIENNA_RNA_PACKAGE_LFOLD_H__


//THIS SECTION WAS MODIFIED BY MAX WARD (2014) AS PART OF HIS HONOURS PROJECT
//It is intended as a drop in replacement for the standard Lfold.h
//Created from a modified version of Vienna 2.14

//A Linked List node
struct mnode {
  char* sstruct;
  float fe;
  int left;
  struct mnode *next;
};

//A modified Lfold function the returns a linked list of local structures, rather than printing them to console
struct mnode* mLfold(const char *string,
            char *structure,
            int maxdist);

//END MODIFICATIONS


/**
 */

/**
 *  \addtogroup local_fold
 *
 *  Local structures can be predicted by a modified version of the
 *  fold() algorithm that restricts the span of all base pairs.
 *  @{
 *    \file Lfold.h
 *    \brief Predicting local MFE structures of large sequences
 *
 *  @}
 */

/**
 *  \addtogroup local_mfe_fold
 *  @{
 *
 *  @}
 */

/**
 *  \brief The local analog to fold().
 * 
 *  Computes the minimum free energy structure including only base pairs
 *  with a span smaller than 'maxdist'
 *
 *  \ingroup local_mfe_fold
 * 
 *  \param string
 *  \param structure
 *  \param maxdist
 */
float Lfold(const char *string,
            char *structure,
            int maxdist);

/**
 *  \brief
 * 
 *  \ingroup local_mfe_fold
 * 
 *  \param string
 *  \param structure
 *  \param maxdist
 *  \param zsc
 *  \param min_z
 */
float Lfoldz( const char *string,
              char *structure,
              int maxdist,
              int zsc,
              double min_z);


/**
 *  \addtogroup local_consensus_fold
 *  @{
 *
 *  @}
 */

/**
 *  \brief
 *
 *  \ingroup local_consensus_fold
 * 
 *  \param strings
 *  \param structure
 *  \param maxdist
 *  \return
 */
float aliLfold( const char **strings,
                char *structure,
                int maxdist);

#endif

import numpy as np
np.set_printoptions(legacy='1.25')
import pandas as pd
import json
import re
import os

from .log_handlers import *

class Yamazaki(object):
    __doc__="""Class for handling data from Tables 1 and 2 in reference article 
    by Yamazaki [1].

    Table 1: Calculated statistical population tensors Bk(J) assuming complete 
    nuclear alignment (Table 1, p.5 [1]).

    Table 2: Calculated angular distribution coefficients for integral spins 
    (Table 2a, p.6-14 [1]) and half-integral spins (Table 2b, p.15-23 [1]).

    References:
        [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    _ROOT = os.path.abspath(os.path.dirname(__file__))
    
    def __init__(self):
        from . import get_data
        data_path = get_data('data')
        data_file = "%s/table1.json"%data_path
        self.data_list = []        
        with open(data_file, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.data_list.append(data_dict)
        jf.close()
        #self.data_list[0]

        # Add Table 2
        yamazaki_table2a = "%s/table2a.json"%data_path
        with open(yamazaki_table2a, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.data_list.append(data_dict)
        jf.close()

        yamazaki_table2b = "%s/table2b.json"%data_path
        with open(yamazaki_table2b, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.data_list.append(data_dict)
        jf.close()
        
    def get_B(self, k, J):
        """Find B for given value of k and J in Table 1 [1].

        Notes:
            Calculated Bk(J) values taken from Table 1 on p.5 of reference 
            article by Yamazaki [1].
        
            [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            k: An integer object representing the order.
            J: A number object (integer or float) representing the spin.

        Returns:
            A float value corresponding to Bk(J) from Table 1 in Yamazaki's 
            paper as described in the `Notes` above.

        Raises:
            Fewer or more than 2 positional arguments ('k', 'J') raises a 
            TypeError exception.

        Example:
            To find B(k=2,J=4):
            > get_B(2,4)
        """
        #print(self.data_list) # instantiation of class creates data list which can be used by all functions in this class using `self.data_list`, `self.data_list[0]`.

        self.k = k
        self.J = J
                
        # Check for k
        if (int(self.k) == 2) or (int(self.k) == 4) or (int(self.k) == 6):

            # Check for J
            if self.J > 0 and self. J < 21:

                try:
                    print("k = {0}".format(self.k))
                    print("J = {0}".format(self.J))

                    MATCH_J = False
                    for jdict in self.data_list[0]:
                        if int(2*float(jdict["J"])) == 2*self.J:
                            MATCH_J = True
                            print("B{0} = {1}".format(int(self.k), jdict["B%i"%int(self.k)]))
                            return float(jdict["B%i"%int(self.k)])

                    if MATCH_J == False:
                        logger.error(f"Only integral or half-integral values of `J` are permitted: J={self.J} is not allowed.")
                        return
                        
                except ValueError:
                    logger.exception("Only integral or half-integral values of J are permitted: \nNumerical values must be given as integers or floats.")
                    return
                except TypeError:
                    logger.exception("Only integral or half-integral values of J are permitted: \nNumerical values must be given as integers or floats.")
                    return
                
            else:
                logger.error("Only the following integral or half-integral values of `J` are permitted: 1<=J<21.")
                if ((2*self.J) %2 == 0) or ((2*self.J) %2 == 1):
                    logger.warning(f"Your value of {self.J} falls outside of this range.")
                else:
                    logger.error(f"J={self.J} is not an acceptable value.")
                return

        else:
            logger.warning(f"Only the following integral values of `k` are permitted: k = 2, 4, or 6.\nYour value of {self.k} falls outside of this range.")
            return

    def get_row_table1(self, J):
        """Method for obtaining individual row data corresponding to the
        calculated statistical population data listed in Table 1 of
        Ref. [1].

        Notes:
            [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            J: A number object (integer or float) representing the spin of the
               nuclear state.

        Returns:
            A list of floating-point objects corresponding to calculated 
            statistical population tensors from Table 1 of Ref. [1].

        Raises:
            Fewer or more than 1 positional arguments ('J') raises a 
            TypeError exception.

        Example:
            To obtain all statistical population tensors for k=2, 4, and 6 
            associated with J=10:

            > get_row_table1(10)
        """
        self.J = J

        if self.J >=1 and self.J < 21:
            if ((2*self.J)%2==0) or ((2*self.J)%2==1):
                for jdict in self.data_list[0]:
                    if int(2*float(jdict["J"])) == int(2*self.J):
                        return [float(jdict["B2"]), float(jdict["B4"]), float(jdict["B6"])]
            else:
                logger.error(f"Only integral or half-integral values of `J` are permitted within range: 1 <= J < 21 \nJ={self.J} is not allowed.")
                return
        else:
            logger.error("Only the following integral or half-integral values of `J` are permitted: 1<=J<21.")
            if ((2*self.J) %2 == 0) or ((2*self.J) %2 == 1):
                logger.warning(f"Your value of {self.J} falls outside of this range.")
            else:
                logger.error(f"J={self.J} is not an acceptable value.")
            return
        
    def get_table1(self,*args):
        """Data table of statistical population tensors for complete alignment:
        B2, B4, and B6 for integral (m=0) spin values 1<=J<=20 and half-integral
        (m=0.5) spin values 3/2<=J<=41/2.

        Notes:
            Calculated Bk(J) values taken from Table 1 on p.5 of reference 
            article by Yamazaki [1].

            [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            None: Both integral (m=0) and half-integral (m=0.5) spin values.
            args: May take 1 additional argument:
                  0: integral (m=0) spin values only;
                  0.5: half-integral (m=0.5) spin values only.

        Returns:
            A DataFrame object containing Bk(J) values from Table 1 in 
            Yamazaki's paper as described in the `Notes` above.  The DataFrame 
            will be populated according to how the function is called (see 
            examples).

        Raises:
            Invalid literals passed to the method raises a ValueError exception.

        Examples:
            (i) Complete table with both m=0 and m=0.5 statistical population 
            tensors:
            > get_table1()

            (ii) Table of m=0 statistical population tensors only:
            > get_table1(0)

            (iii) Table of m=0.5 statistical population tensors only:
            > get_table1(0.5)
        """
        stat_tensor_list = self.data_list[0]
        stat_tensor_df = pd.DataFrame(stat_tensor_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        if len(args) == 0:
            return stat_tensor_df
        elif len(args) == 1:
            try:
                m = args[0]
                if int(2*m) == 0:
                    #return stat_tensor_df.loc[stat_tensor_df['m']=='0']
                    return stat_tensor_df.loc[stat_tensor_df['m']==0]
                elif int(2*m) == 1:
                    #return stat_tensor_df.loc[stat_tensor_df['m']=='0.5']
                    return stat_tensor_df.loc[stat_tensor_df['m']==0.5]
                else:
                    logger.error("Argument not accepted.  Call the function using: \n get_B_table(0) # m=0 results only \n get_B_table(0.5) # m=0.5 results only")
            except ValueError:
                logger.exception("Incorrect input argument.  Call the function using: \n get_B_table(0) # m=0 results only \n get_B_table(0.5) # m=0.5 results only") 
            return
        else:
            logger.error("Wrong input arguments!  Call the function using: \n get_B_table() # m=0 and m=0.5 results \n get_B_table(0) # m=0 results only \n get_B_table(0.5) # m=0.5 results only")
            return

    def get_table2a(self):
        """Table of angular distribution coefficients associated with 
        even-integral spins given in Table 2(a) of Ref. [1].

        Notes:
            [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding angular
            distribution coefficients Fk, FkBk, and Uk listed in
            Table 2(a) [1] for even-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_table2a()
        """
        table2a_list = self.data_list[1]
        table2a_df = pd.DataFrame(table2a_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return table2a_df

    def get_table2b(self):
        """Table of angular distribution coefficients associated with 
        odd-integral spins given in Table 2(b) of Ref. [1].

        Notes:
            [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding angular
            distribution coefficients Fk, FkBk, and Uk listed in
            Table 2(b) [1] for odd-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_table2b()
        """
        table2b_list = self.data_list[2]
        table2b_df = pd.DataFrame(table2b_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return table2b_df
    
    def get_row_table2(self, Ji, Jf, L1, L2, *args, **kwargs):
        """Method for manipulating individual row data corresponding to the
        angular distribution coefficients listed in Tables 2(a) and (b) of
        Ref. [1].

        Notes:
            [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            Ji: A number object (integer or float) representing the initial 
                spin.
            Jf: A number object (integer or float) representing the final spin.
            L1: An integer object representing the first multipole order.
            L2: An integer object representing the second multipole order.

            args: Optional argument to specify the k-order.
                  Only integers of 2 or 4 are acceptable values for k.
            kwargs: Optional keyword argument to specify angular distribution
                    coefficient as:
                    'coeff'='F'
                    'coeff'='BF'
                    'coeff'='U'

        Returns:
            A list of floating-point objects or a single float corresponding to
            angular distribution coefficients from Tables 2(a) and (b) of Ref.
            [1967Ya05] depending on the input arguments passed to the method.

        Raises:
            Invalid literals passed as arguments to the method raises a 
            ValueError exception.

            Incorrect data types passed as arguments to the method raises a 
            TypeError exception.

        Examples:
            For the transition associated with Ji=15.5, Jf=18.5, L1=3, L2=3:-

            (i) To get a complete list all coefficients in corresponding row:
            > get_row_table2(15.5,18.5,3,3)

            (ii) To get a list of all k=2 coefficients in row:
            > get_row_table2(15.5,18.5,3,3,2)

            (iii) To get a list of all 'F' coefficients in row:
            > get_row_table2(15.5,18.5,3,3,coeff='F')

            (iv) To get the 'U4' coefficient from the row:
            > get_row_table2(15.5,18.5,3,3,4,coeff='U')
        """
        self.Ji = Ji
        self.Jf = Jf
        self.L1 = L1
        self.L2 = L2
        self.args = args
        self.kwargs = kwargs

        k = None
        if len(args) == 0:
            k = None
        elif len(args) == 1:
            try:
                k = int(args[0])
            except ValueError:
                logger.exception(f"Wrong value!  Argument must be a bytes-like object or a number; {args[0]} is not an acceptable value.  Only integral 'k' values of 2 or 4 are acceptable.")
                return
            except TypeError:
                logger.exception(f"Wrong type!  Argument must be a bytes-like object or a number, not a {type(args[0])}.  Only integral 'k' values of 2 or 4 are acceptable.")
                return

        C = None
        if len(kwargs) == 0 or kwargs == {}:
            C = None
        else:
            for key, value in kwargs.items():
                if key == "coeff":
                    C = value
                    if type(C) is str:
                        C = C.upper()
                else:
                    logger.error("Key not recognized: Use 'coeff' as key argument.")
                    return

        # Check for even and half-integral J and ranges on Ji and Jf
        RANGE_VALID_Ji = False
        RANGE_VALID_Jf = False
        table2_data = None
        if ((2*self.Ji)%2 == 0) and ((2*self.Jf)%2 == 0):
            table2_data = self.data_list[1]
            if self.Ji >= 1 and self.Ji <= 15:
                RANGE_VALID_Ji = True
            if self.Jf >= 0 and self.Jf <= 18:
                RANGE_VALID_Jf = True
            
        elif ((2*self.Ji)%2 == 1) and ((2*self.Jf)%2 == 1):
            table2_data = self.data_list[2]
            if self.Ji >= 1.5 and self.Ji <= 15.5:
                RANGE_VALID_Ji = True
            if self.Jf >= 0.5 and self.Jf <= 18.5:
                RANGE_VALID_Jf = True

        if RANGE_VALID_Ji == False or RANGE_VALID_Jf == False:
            if RANGE_VALID_Ji == False:
                logger.error(f"'Ji'={self.Ji} outside range of permissive values:\n 1 <= Ji <= 15 (integral J) \n 3/2 <= Ji <= 31/2 (half-integral J)")
            if RANGE_VALID_Jf == False:
                logger.error(f"'Jf'={self.Jf} outside range of permissive values:\n 0 <= Jf <= 18 (integral J) \n 1/2 <= Jf <= 37/2 (half-integral J)")
            return

        # Check L-values
        L1_VALID = False
        if self.L1 > 0 and self.L1 <= 3:
            L1_VALID = True
        else:
            logger.warning(f"L1={self.L1} is not a permissive multipole order.")
            
        L2_VALID = False
        if self.L2 > 0 and self.L2 <= 3:
            L2_VALID = True
        else:
            logger.warning(f"L2={self.L2} is not a permissive multipole order.")

        if L1_VALID == False or L2_VALID == False:
            logger.error("Multipole order must fall within range: 0 < L <= 3")
            return
                        
        for jdict in table2_data:
            if int(2*float(jdict["Ji"])) == int(2*self.Ji) and int(2*float(jdict["Jf"])) == int(2*self.Jf) and int(2*float(jdict["L1"])) == int(2*self.L1) and int(2*float(jdict["L2"])) == int(2*self.L2) and RANGE_VALID_Ji == True and RANGE_VALID_Jf == True and L1_VALID == True and L2_VALID == True:

                WRONG_k = False
                if k is not None:
                    if k % 2 == 1:
                        WRONG_k = True
                    if (k < 2) or (k > 4):
                        WRONG_k = True
                
                if k == None and C == None:
                    # Returns complete row
                    return [jdict["F2"], jdict["B2F2"], jdict["F4"], jdict["B4F4"], jdict["U2"], jdict["U4"]]
                elif (k == 2 or k == 4) and C == None:
                    # Returns all specified k values in row
                    return [jdict["F%i"%k], jdict["B%iF%i"%(k,k)], jdict["U%i"%k]]
                elif k == None and (C == "F" or C == "BF" or C == "U"):
                    # Returns all specified coefficients in row
                    if len(C) == 1:
                        return [jdict["%s2"%C], jdict["%s4"%C]]
                    elif len(C) == 2:
                        return [jdict["B2F2"], jdict["B4F4"]]

                elif (k == 2 or k == 4) and (C == "F" or C == "BF" or C == "U"):
                    # Returns all specified coefficients in row of k value
                    if len(C) == 1:
                        return jdict["%s%i"%(C,k)]
                    elif len(C) == 2:
                        return jdict["B%iF%i"%(k,k)]
                else:
                    if WRONG_k == True:
                        logger.error("Argument `k` must be given as a value of 2 or 4 only.")
                    else:
                        logger.error("Keyword argument `coeff` only takes string values of 'F', 'BF', or 'U'.")
                    return


    def yamazaki2file(self,table,format):
        """Convert and dump data from Tables 1, 2(a), and 2(b) of Ref. [1] into 
        a CSV or JSON formatted file.

        Notes:
            [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            table: String argument representing desired table printout:-
                   'T1': Table 1 [1]
                   'T2A': Table 2(a) [1]
                   'T2B': Table 2(b) [1]
            format: String argument to indicate preferred file format:-
                   'CSV': Comma Separated Value format.
                   'JSON': JavaScript Object Notation format.

        Returns:
            None.

        Examples:
            To print Table 1 to file in CSV format:
            > yamazaki2file('T1','CSV') # dumps `YamazakiTable1.csv` in pwd

            To print Table 2(a) to file in JSON format:
            > yamazaki2file('T2A','JSON') # dumps `YamazakiTable2a.json` in pwd
        """
        self.table = table
        self.format = format

        tables = {"T1":[self.data_list[0], "YamazakiTable1"],
                  "T2A":[self.data_list[1], "YamazakiTable2a"],
                  "T2B":[self.data_list[2], "YamazakiTable2b"]}
        formats = ["CSV","JSON"]

        FILE_FORMAT = False
        TABLE_YAMAZAKI = False
        PRINT_PROBLEM = False
        
        for f in formats:
            if f == self.format:
                FILE_FORMAT = True
                for key, value in tables.items():
                    if key == self.table:
                        TABLE_YAMAZAKI = True
                        table_data = value[0]
                        with open("%s.%s"%(value[1],f.lower()), mode="w") as outfile:
                            if f.upper() == "JSON":
                                json.dump(table_data, outfile, indent=4, ensure_ascii=False)
                                logger.info("{0}.{1} printed to file in {2}".format(value[1],f.lower(),os.getcwd()))
                                outfile.close()
                            elif f.upper() == "CSV":
                                table_df = pd.DataFrame(table_data)
                                table_df.to_csv(outfile, index=False)
                                logger.info("{0}.{1} printed to file in {2}".format(value[1],f.lower(),os.getcwd()))
                                outfile.close()

        if FILE_FORMAT == False:
            PRINT_PROBLEM = True
            logger.error("File format not handled: Specify 'CSV' or 'JSON'")
            for key in tables.keys():
                if key == self.table: TABLE_YAMAZAKI = True
        if TABLE_YAMAZAKI == False:
            PRINT_PROBLEM = True
            logger.error("Table from Yamazaki's paper not correctly specified: \n 'T1' - Table 1; 'T2A' - Table 2(a); 'T2B' - Table 2(b).")

        if PRINT_PROBLEM == True:
            logger.error("File not printed.")


class RoseAndBrink(Yamazaki):
    __doc__="""Class for handling tabulated data in Appendex of the review 
    article by Rose and Brink [2].

    (i)   Angular distribution Rk(L1 L2 Ji Jf) coefficients (L2=L1+1) for 
          integral and half-integral spins.

    (ii)  Angular distribution Uk(L1 Ji Jf) and Uk(L2 Ji Jf) coefficients 
          (L2=L1+1) for integral and half-integral spins.
    
    (iii) Population tensor Sk(l1 l2 J s) coefficients for integral and 
          half-integral spins.
    
    (iv)  Statistical tensor p(J m) coefficients for integral and half-integral
          spins.

    References:
        [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
            p. 306 (1967).
    """
    _ROOT = os.path.abspath(os.path.dirname(__file__))
    
    def __init__(self):
        from . import get_data
        data_path = get_data('data')
        self.rose_brink_list = []

        # [0]: R (integral), [1]: R (half-integral)
        # [2]: U (integral and half-integral together)
        # [3]: S (integral), [4]: S (half-integral)
        # [5]: p (integral), [6]: p (half-integral)
        
        # Handle Rose and Brink data tables
        # R-coefficients
        rose_brink_Ra = "%s/rose_brink_coeff_R_a.json"%data_path
        with open(rose_brink_Ra, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.rose_brink_list.append(data_dict)
        jf.close()

        rose_brink_Rb = "%s/rose_brink_coeff_R_b.json"%data_path
        with open(rose_brink_Rb, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.rose_brink_list.append(data_dict)
        jf.close()

        # U-coefficients
        rose_brink_U = "%s/rose_brink_coeff_U.json"%data_path
        with open(rose_brink_U, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.rose_brink_list.append(data_dict)
        jf.close()

        # S-coefficients
        rose_brink_Sa = "%s/rose_brink_coeff_S_a.json"%data_path
        with open(rose_brink_Sa, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.rose_brink_list.append(data_dict)
        jf.close()

        rose_brink_Sb = "%s/rose_brink_coeff_S_b.json"%data_path
        with open(rose_brink_Sb, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.rose_brink_list.append(data_dict)
        jf.close()

        # p-coefficients
        rose_brink_pa = "%s/rose_brink_coeff_p_a.json"%data_path
        with open(rose_brink_pa, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.rose_brink_list.append(data_dict)
        jf.close()

        rose_brink_pb = "%s/rose_brink_coeff_p_b.json"%data_path
        with open(rose_brink_pb, mode='r') as jf:
            data_dict = json.loads(jf.read())
            self.rose_brink_list.append(data_dict)
        jf.close()

    def get_tableRa(self):
        """Table of angular distribution coefficients associated with 
        even-integral spins for the Rk(L1 L2 Ji Jf) coefficients, where L2=L1+1,
        given in the Appendix of the Rose and Brink review article [2].

        Notes:
            [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
                p. 306 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding Rk angular
            distribution coefficients for k=2,4,6, and 8 as tabulated in the 
            Appendix of Ref. [2] for even-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_tableRa()
        """
        tableRa_list = self.rose_brink_list[0]
        tableRa_df = pd.DataFrame(tableRa_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return tableRa_df

    def get_tableRb(self):
        """Table of angular distribution coefficients associated with 
        half-integral spins for the Rk(L1 L2 Ji Jf) coefficients, where L2=L1+1,
        given in the Appendix of the Rose and Brink review article [2].

        Notes:
            [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
                p. 306 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding Rk angular
            distribution coefficients for k=2,4,6, and 8 as tabulated in the 
            Appendix of Ref. [2] for half-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_tableRb()
        """
        tableRb_list = self.rose_brink_list[1]
        tableRb_df = pd.DataFrame(tableRb_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return tableRb_df

    def get_tableU(self):
        """Table of angular distribution coefficients associated with both  
        integral and half-integral spins for the Uk(L1 Ji Jf) and Uk(L2 Ji Jf) 
        coefficients, where L2 = L1 + 1, given in the Appendix of the 
        Rose and Brink review article [2].

        Notes:
            [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
                p. 306 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding Uk angular
            distribution coefficients for k=2,4,6, and 8 as tabulated in the 
            Appendix of Ref. [2] for integral and half-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_tableRb()
        """
        tableU_list = self.rose_brink_list[2]
        tableU_df = pd.DataFrame(tableU_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return tableU_df

    def get_tableSa(self):
        """Table of population-tensor coefficients associated with 
        even-integral spins for the Sk(l1 l2 J s) coefficients as given in the 
        Appendix of the Rose and Brink review article [2].

        Notes:
            [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
                p. 306 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding Sk population 
            tensor coefficients for k=0,2,4,6, and 8 as tabulated in the 
            Appendix of Ref. [2] for even-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_tableSa()
        """
        tableSa_list = self.rose_brink_list[3]
        tableSa_df = pd.DataFrame(tableSa_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return tableSa_df

    def get_tableSb(self):
        """Table of population-tensor coefficients associated with 
        half-integral spins for the Sk(l1 l2 J s) coefficients as given in the 
        Appendix of the Rose and Brink review article [2].

        Notes:
            [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
                p. 306 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding Sk population 
            tensor coefficients for k=0,2,4,6, and 8 as tabulated in the 
            Appendix of Ref. [2] for half-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_tableSb()
        """
        tableSb_list = self.rose_brink_list[4]
        tableSb_df = pd.DataFrame(tableSb_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return tableSb_df

    def get_tablePa(self):
        """Table of statistical-tensor coefficients associated with 
        even-integral spins for the pk(J m) coefficients as given in the 
        Appendix of the Rose and Brink review article [2].

        Notes:
            [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
                p. 306 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding pk statistical 
            tensor coefficients for even-k<=8 and magnetic substate projections 
            m=0,1,2,3,4,5, and 6 as tabulated in the Appendix of Ref. [2] for 
            even-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_tablePa()
        """
        tablePa_list = self.rose_brink_list[5]
        tablePa_df = pd.DataFrame(tablePa_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return tablePa_df

    def get_tablePb(self):
        """Table of statistical-tensor coefficients associated with 
        half-integral spins for the pk(J m) coefficients as given in the 
        Appendix of the Rose and Brink review article [2].

        Notes:
            [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
                p. 306 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding pk statistical 
            tensor coefficients for even-k<=8 and magnetic substate projections 
            m=1/2,3/2,5/2,7/2,9/2,11/2, and 13/2 as tabulated in the Appendix 
            of Ref. [2] for half-integral spins.

        Raises:
            Passing arguments to the function raises a TypeError exception.

        Example:
            > get_tablePb()
        """
        tablePa_list = self.rose_brink_list[6]
        tablePa_df = pd.DataFrame(tablePa_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return tablePa_df

    def rosebrink2file(self,table,format):
        """Convert and dump data from tables in the Appendix of the Rose and 
        Brink review article [2] into a CSV or JSON formatted file.

        Notes:
            [2] H.J. Rose and D.M. Brink, Rev. Mod. Phys., Vol. 39, Num. 2, 
                p. 306 (1967).

        Arguments:
            table: String argument representing desired table printout:-
                   'RA': Rk(L1 L2 Ji Jf) integral-J [2]
                   'RB': Rk(L1 L2 Ji Jf) half-integral-J [2]
                   'U' : Uk(L1 Ji Jf) and Uk(L2 Ji Jf) [2]
                   'SA': Sk(l1 l2 J s) integral-J [2]
                   'SB': Sk(l1 l2 J s) half-integral-J [2]
                   'PA': pk(J m) integral-J [2]
                   'PB': pk(J m) half-integral-J [2]
            format: String argument to indicate preferred file format:-
                   'CSV': Comma Separated Value format.
                   'JSON': JavaScript Object Notation format.

        Returns:
            None.

        Examples:
            To print Table 1 to file in CSV format:
            > rosebrink2file('RA','CSV') # dumps `RoseBrinkTableRa.csv` in pwd

            To print Table 2(a) to file in JSON format:
            > rosebrink2file('U','JSON') # dumps `RoseBrinkTableU.json` in pwd
        """
        self.table = table
        self.format = format

        tables = {"RA":[self.rose_brink_list[0], "RoseBrinkTableRa"],
                  "RB":[self.rose_brink_list[1], "RoseBrinkTableRb"],
                  "U":[self.rose_brink_list[2], "RoseBrinkTableU"],
                  "SA":[self.rose_brink_list[3], "RoseBrinkTableSa"],
                  "SB":[self.rose_brink_list[4], "RoseBrinkTableSb"],
                  "PA":[self.rose_brink_list[5], "RoseBrinkTablePa"],
                  "PB":[self.rose_brink_list[6], "RoseBrinkTablePb"]}
        formats = ["CSV","JSON"]

        FILE_FORMAT = False
        TABLE_ROSE_BRINK = False
        PRINT_PROBLEM = False
        
        for f in formats:
            if f == self.format:
                FILE_FORMAT = True
                for key, value in tables.items():
                    if key == self.table:
                        TABLE_ROSE_BRINK = True
                        table_data = value[0]
                        with open("%s.%s"%(value[1],f.lower()), mode="w") as outfile:
                            if f.upper() == "JSON":
                                json.dump(table_data, outfile, indent=4, ensure_ascii=False)
                                logger.info("{0}.{1} printed to file in {2}".format(value[1],f.lower(),os.getcwd()))
                                outfile.close()
                            elif f.upper() == "CSV":
                                table_df = pd.DataFrame(table_data)
                                table_df.to_csv(outfile, index=False)
                                logger.info("{0}.{1} printed to file in {2}".format(value[1],f.lower(),os.getcwd()))
                                outfile.close()

        if FILE_FORMAT == False:
            PRINT_PROBLEM = True
            logger.error("File format not handled: Specify 'CSV' or 'JSON'")
            for key in tables.keys():
                if key == self.table: TABLE_ROSE_BRINK = True
        if TABLE_ROSE_BRINK == False:
            PRINT_PROBLEM = True
            logger.error("Table from Rose and Brink's paper not correctly specified: \n 'RA' - Table R (integral J); 'RB' - Table R (half-integral J); \n 'U' - Table U; \n 'SA' - Table S (integral J); 'SB' - Table S (half-integral J); \n 'PA' - Table p (integral J); 'PB' - Table p (half-integral J).")

        if PRINT_PROBLEM == True:
            logger.error("File not printed.")    

import numpy as np
import pandas as pd
import json
import re
import os

class Tables(object):
    __doc__="""Class for handling data from Tables 1 and 2 in Ref. [1967Ya05].

    Table 1: Calculated statistical population tensors Bk(J) assuming complete 
    nuclear alignment (Table 1, p.5 [1967Ya05]).

    Table 2: Calculated angular distribution coefficients for integral spins 
    (Table 2a, p.6-14 [1967Ya05]) and half-integral spins (Table 2b, p.15-23 
    [1967Ya05]).

    References:
    [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
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
        """Find B for given value of k and J in Table 1 [1967Ya05].

        Notes:
            Calculated Bk(J) values taken from Table 1 on p.5 of reference 
            article [1967Ya05].
        
            [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            k: An integer object representing the order.
            J: A number object (integer or float) representing the spin.

        Returns:
            A float value corresponding to Bk(J) from Table 1 in Yamazaki's 
            paper as described in the `Notes` above.

        Example:
            To find B(k=2,J=4):
            > get_B(2,4)
        """
        #print(self.data_list) # instantiation of class creates data list which can be used by all functions in this class using `self.data_list`, `self.data_list[0]`.

        self.k = k
        self.J = J
                
        # Check for k
        if (int(self.k == 2)) or (int(self.k == 4)) or (int(self.k == 6)):

            # Check for J
            if self.J > 0 and self. J < 21:

                # Check for integral values of J first
                try:
                    print("k = {0}".format(self.k))
                    print("J = {0}".format(self.J))

                    for jdict in self.data_list[0]:
                        if int(2*float(jdict["J"])) == int(2*self.J):
                            print("B{0} = {1}".format(int(self.k), jdict["B%i"%int(self.k)]))
                            return float(jdict["B%i"%int(self.k)])
                        
                except ValueError:
                    print("Only integral or half-integral values of J are permitted.")
                    print("Not J=",self.J)
                    return
                
            else:
                print("Only the following integral values of `J` are permitted:")
                print("1<=J<21.")
                print("Your value of {0} does not fall within this range.".format(self.J))
                return

        else:
            print("Only the following integral values of `k` are permitted:")
            print("k = 2, 4, or 6.")
            print("Your value of {0} does not satisfy this requirement.".format(self.k))
            return

    def get_row_table1(self, J):
        """Method for obtaining individual row data corresponding to the
        calculated statistical population data listed in Table 1 of
        Ref. [1967Ya05].

        Notes:
            [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            J: A number object (integer or float) representing the spin of the
               nuclear state.

        Returns:
            A list of floating-point objects corresponding to calculated 
            statistical population tensors from Table 1 of Ref. [1967Ya05].

        Example:
            To obtain all statistical population tensors for k=2,4, and 6 
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
                print("Spin must be integral of half integral within range: 1 <= J < 21")
                return
        else:
            print("Spin must be integral of half integral within range: 1 <= J < 21")
            return
        
    def get_table1(self,*args):
        """Data table of statistical population tensors for complete alignment:
        B2, B4, and B6 for integral (m=0) spin values 1<=J<=20 and half-integral
        (m=0.5) spin values 3/2<=J<=41/2.

        Notes:
            Calculated Bk(J) values taken from Table 1 on p.5 of reference 
            article by T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            None: Both integral (m=0) and half-integral (m=0.5) spin values.
            args: May take 1 additional argument:
                  0: integral (m=0) spin values only;
                  0.5: half-integral (m=0.5) spin values only.

        Returns:
            A DataFrame object containing Bk(J) values from Table 1 in 
            Yamazaki's paper as described in the `Notes` above.  The DataFrame 
            will be populated according to how the function is called:

            (i)   Both m=0 and m=0.5 statistical population tensors:
                  get_B_table()
            (ii)  m=0 statistical population tensors:
                  get_B_table(0)
            (iii) m=0.5 statistical population tensors:
                  get_B_table(0.5)
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
                    print("Argument not accepted.")
                    print("Please the function using:")
                    print(" get_B_table(0) # m=0 results only")
                    print(" get_B_table(0.5) # m=0.5 results only")
            except ValueError:
                print("Incorrect input argument:")
                print("Call the function using:")
                print(" get_B_table(0) # m=0 results only")
                print(" get_B_table(0.5) # m=0.5 results only") 
            return
        else:
            print("Wrong input arguments!")
            print("Call the function using:")
            print(" get_B_table() # m=0 and m=0.5 results")
            print(" get_B_table(0) # m=0 results only")
            print(" get_B_table(0.5) # m=0.5 results only")
            return

    def get_table2a(self):
        """Table of angular distribution coefficients associated with even-integral
        spins given in Table 2(a) of Ref. [1967Ya05].

        Notes:
            [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding angular
            distribution coefficients Fk, FkBk, and Uk listed in
            Table 2(a) [1967Ya05] for even-integral spins.

        Usage:
            > get_table2a()
        """
        table2a_list = self.data_list[1]
        table2a_df = pd.DataFrame(table2a_list)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_row', None)

        return table2a_df

    def get_table2b(self):
        """Table of angular distribution coefficients associated with odd-integral
        spins given in Table 2(b) of Ref. [1967Ya05].

        Notes:
            [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            None.

        Returns:
            A DataFrame object containing the corresponding angular
            distribution coefficients Fk, FkBk, and Uk listed in
            Table 2(b) [1967Ya05] for odd-integral spins.

        Usage:
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
        Ref. [1967Ya05].

        Notes:
            [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            Ji: A number object (integer or float) representing the initial spin.
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
                print("Argument must be integral:") 
                print("Only 'k' values of 2 or 4 are acceptable.")
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
                    print("Key not recognized: Use 'coeff' as key argument.")

        table2_data = None
        if ((2*self.Ji)%2 == 0) and ((2*self.Jf)%2 == 0):
            table2_data = self.data_list[1]
        elif ((2*self.Ji)%2 == 1) and ((2*self.Jf)%2 == 1):
            table2_data = self.data_list[2]
            
        for jdict in table2_data:
            if int(2*float(jdict["Ji"])) == int(2*self.Ji) and int(2*float(jdict["Jf"])) == int(2*self.Jf) and int(2*float(jdict["L1"])) == int(2*self.L1) and int(2*float(jdict["L2"])) == int(2*self.L2):
                    
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
                    print("Wrong argument value and/or keyword argument:")
                    print("Use integral values of 2 or 4 only for k args")
                    print("Keyword argument `coeff` only takes string values of 'F', 'BF', or 'U'")
                    return


    def table2file(self,table,format):
        """Convert and dump data from Tables 1, 2(a), and 2(b) of Ref. [1967Ya05]
        into a CSV or JSON formatted file.

        Notes:
            [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            table: String argument representing desired table printout:-
                   'T1': Table 1 [1967Ya05]
                   'T2A': Table 2(a) [1967Ya05]
                   'T2B': Table 2(b) [1967Ya05]
            format: String argument to indicate preferred file format:-
                   'CSV': Comma Separated Value format.
                   'JSON': JavaScript Object Notation format.

        Returns:
            None.

        Examples:
            To print Table 1 to file in CSV format:
            > table2file('T1','CSV') # dumps `yamazaki_table1.csv` in pwd

            To print Table 2(a) to file in JSON format:
            > table2file('T2A','JSON') # dumps `yamazaki_table2a.json` in pwd
        """
        self.table = table
        self.format = format

        tables = {"T1":[self.data_list[0], "yamazaki_table1"], "T2A":[self.data_list[1], "yamazaki_table2a"],
                  "T2B":[self.data_list[2], "yamazaki_table2b"]}
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
                                outfile.close()
                            elif f.upper() == "CSV":
                                table_df = pd.DataFrame(table_data)
                                table_df.to_csv(outfile, index=False)

        if FILE_FORMAT == False:
            PRINT_PROBLEM = True
            print("File format not handled: Specify 'CSV' or 'JSON'")
            for key in tables.keys():
                if key == self.table: TABLE_YAMAZAKI = True
        if TABLE_YAMAZAKI == False:
            PRINT_PROBLEM = True
            print("Table from Yamazaki's paper not correctly specified:")
            print("'T1' - Table 1; 'T2A' - Table 2(a); 'T2B' - Table 2(b).")

        if PRINT_PROBLEM == True:
            print("File not printed.")
        
        

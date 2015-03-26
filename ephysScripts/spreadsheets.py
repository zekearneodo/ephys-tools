__author__ = 'zeke'
'''
Class for handling google data spreadsheets.

'''

import gdata.spreadsheet.service as gss
import gdata


class Spreadsheets:
    # a collection of spreadsheets
    def __init__(self, client):
        # The client is given and is logged in
        self.client = client
        self.spreadsheet_key = ''
        self.worksheet_id = ''
        self.cells_feed = ''
        pass

    def list_spreadsheets(self, verbose=False):
        # gets the feed of worksheets and prints the list
        """
        returns a dictionary containing keys for titles and id of the spreadshet as value:
        (it helps when you want to load a particular sheet)
        sheets[title]=sheetId
        :param verbose: if True, it will print the list of titles -> id.
                default is False
        """
        feed = self.client.GetSpreadsheetsFeed()
        sheets = {}
        for i, entry in enumerate(feed.entry):
            sheets[entry.title.text] = entry.id.text.rsplit('/', 1)[1]
            if verbose:
                print '{0:40} -> {1}'.format(entry.title.text, sheets[entry.title.text])

        return sheets

    def _list_worksheets(self, verbose=False):
        # gets the feed of worksheets and prints the list
        """
        returns a dictionary containing keys for titles and id of the spreadshet as value:
        (it helps when you want to load a particular sheet)
        sheets[title]=sheetId
        :param verbose: if True, it will print the list of titles -> id.
                default is False
        """
        if not self.spreadsheet_key:
            print "No Spreadsheet selected"
            return
        else:
            spreadsheet = self.client.GetWorksheetsFeed(self.spreadsheet_key)
            sheets = {}
            for i, entry in enumerate(spreadsheet.entry):
                sheets[entry.title.text] = entry.id.text.rsplit('/', 1)[1]
                if verbose:
                    print '{0:40} -> {1}'.format(entry.title.text, sheets[entry.title.text])
            return sheets

    def get_spreadsheet(self, title):
        """
        returns the sheet with the entered title
        sets it as the current spreadsheet in self.spreadsheet_key
        """
        sheets = self.list_spreadsheets(verbose=False)

        try:
            self.spreadsheet_key = sheets[title]
            sheet = self.client.GetWorksheetsFeed(self.spreadsheet_key)
            return sheet
        except:
            print 'Spreadsheet {0} not found'.format(title)
            sheet_id = ''
            sheet = {}

    def get_worksheet(self, title=''):
        """
        gets a worksheet (a cells_feed) from the spreadsheet that is pointed at in self._spreadsheet_key
        sets it
        :param title: (optional) the title of the worksheet. The default is to use the one that is in self.worksheet_id
        """
        sheets = self._list_worksheets()
        if not title:
            if not self.worksheet_id:
                print 'No title entered and no title of worksheet previously set'
                return None
        else:
            self.worksheet_id = sheets[title]

        try:
            worksheet = self.client.GetCellsFeed(self.spreadsheet_key, self.worksheet_id)
            self.cells_feed = worksheet
            return worksheet
        except:
            print 'Worksheet {0} not found'.format(self.worksheet_id)
            return None

    def get_col_values(self, col, row_range=[]):

        """
        Get a list with the text values of a column in a worksheet

        :param col: column number (int)
        :param row_range: range of rows, starting from 1.
                e.g: [2,5] gives from 2 to 5 includingly
        :return: col_content
                list of the values.
                Empty cells have NoneType text, which is translated to empty string.
        """
        if not self.cells_feed:
            print "No cells feed loaded"
            return None

        if not row_range:
            row_range = [1, 'end']

        # make a query for the column
        if row_range[1] is 'end':
            row_range[1] = int(self.cells_feed.row_count.text)

        query = gss.CellQuery()
        query['min-col'] = str(col)
        query['max-col'] = str(col)
        query['min-row'] = str(row_range[0])
        query['max-row'] = str(row_range[1])
        query['return-empty'] = 'true'
        col_feed = self.client.GetCellsFeed(self.spreadsheet_key, self.worksheet_id, query=query)
        col_content = []
        for j, row in enumerate(col_feed.entry):
            content = row.content.text
            # if content is None:
            #     content = ''
            col_content.append(content)

        return col_content

    def get_row_values(self, row, col_range=[]):

        """
        Get a list with the text values of a column in a worksheet

        :param row: row number (int)
        :param col_range: range of cols, starting from 1.
                e.g: [2,5] gives from 2 to 5 includingly
        :return: row_content
                list of the values.
                Empty cells have NoneType text, which is translated to empty string.
        """
        if not self.cells_feed:
            print "No cells feed loaded"
            return None

        if not col_range:
            col_range = [1, 'end']

        # make a query for the column
        if col_range[1] is 'end':
            col_range[1] = int(self.cells_feed.col_count.text)

        #make a query for the column
        query = gss.CellQuery()
        query['min-row'] = str(row)
        query['max-row'] = str(row)
        query['min-col'] = str(col_range[0])
        query['max-col'] = str(col_range[1])
        query['return-empty'] = 'true'
        row_feed = self.client.GetCellsFeed(self.spreadsheet_key, self.worksheet_id, query=query)
        row_content = []
        for j, row in enumerate(row_feed.entry):
            content = row.content.text
            # if content is None:
            #     content = ''
            row_content.append(content)

        return row_content

    def get_cell_value(self, row, col):
        """
        make a query for the content of a cell
        :param row:
        :param col:
        :return:
        """
        if not self.cells_feed:
            print "No cells feed loaded"
            return None

        # make a query for the cell
        cell_address = 'R{0}C{1}'.format(row, col)
        cell_feed = self.client.GetCellsFeed(self.spreadsheet_key, self.worksheet_id, cell=cell_address)
        content = cell_feed.content.text
        # if content is None:
        #     content = ''

        return content

    def get_count(self):
        if not self.cells_feed:
            print "No cells feed loaded"
            return None

        return [int(self.cells_feed.row_count.text), int(self.cells_feed.col_count.text)]

    def set_cell_value(self, row, col, value):
        if not self.spreadsheet_key or not self.worksheet_id:
            print "No spreadsheet or worksheet specified"
            return None

        entry = self.client.UpdateCell(row=row, col=col, inputValue=value,
                                          key=self.spreadsheet_key, wksht_id=self.worksheet_id)
        if isinstance(entry, gdata.spreadsheet.SpreadsheetsCell):
            return 0
        else:
            print "Error writing to cell; cell not written"
            return 1



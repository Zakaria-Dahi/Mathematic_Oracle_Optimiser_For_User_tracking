function [status, message]=xlwrite(filename,A,sheet, range)

if exist('org.apache.poi.ss.usermodel.WorkbookFactory', 'class') ~= 8 ...
    || exist('org.apache.poi.hssf.usermodel.HSSFWorkbook', 'class') ~= 8 ...
    || exist('org.apache.poi.xssf.usermodel.XSSFWorkbook', 'class') ~= 8
    
    error('xlWrite:poiLibsNotLoaded',...
        'The POI library is not loaded in Matlab.\nCheck that POI jar files are in Matlab Java path!');
end

import org.apache.poi.ss.usermodel.*;
import org.apache.poi.hssf.usermodel.*;
import org.apache.poi.xssf.usermodel.*;
import org.apache.poi.ss.usermodel.*;

import org.apache.poi.ss.util.*;

status=0;

 
if nargin < 3; sheet = []; end
if nargin < 4; range = []; end

 
if nargin < 4 && ~isempty(strfind(sheet,':'))
    range = sheet;
    sheet = [];
end

 
if isempty(A)
    error('xlwrite:EmptyInput', 'Input array is empty!');
end
 
if ndims(A) > 2
	error('xlwrite:InputDimension', ...
        'Dimension of input array should not be higher than two.');
end

 
java.lang.System.setProperty('user.dir', pwd);

 
xlsFile = java.io.File(filename);


if xlsFile.isFile()
    fileIn = java.io.FileInputStream(xlsFile);
    xlsWorkbook = WorkbookFactory.create(fileIn);
else
    [~,~,fileExt] = fileparts(filename);
    
    switch lower(fileExt)
        case '.xls'
            xlsWorkbook = HSSFWorkbook();
        case '.xlsx'
            xlsWorkbook = XSSFWorkbook();
        otherwise
            xlsWorkbook = XSSFWorkbook();
            
            filename = [filename '.xlsx'];
    end
end


if ~isempty(sheet)
    if isnumeric(sheet)
        if xlsWorkbook.getNumberOfSheets() >= sheet && sheet >= 1
            xlsSheet = xlsWorkbook.getSheetAt(sheet-1);
        else
            xlsSheet = [];
        end
    else
        xlsSheet = xlsWorkbook.getSheet(sheet);
    end
    
    if isempty(xlsSheet)
        warning('xlwrite:AddSheet', 'Added specified worksheet.');
        
        if isnumeric(sheet)
            xlsSheet = xlsWorkbook.createSheet(['Sheet ' num2str(sheet)]);
        else
            sheet = WorkbookUtil.createSafeSheetName(sheet);
            xlsSheet = xlsWorkbook.createSheet(sheet);
        end
    end
    
else
    nSheets = xlsWorkbook.getNumberOfSheets();
    
    if nSheets < 1
        xlsSheet = xlsWorkbook.createSheet('Sheet 1');
    else
        xlsSheet = xlsWorkbook.getSheetAt(0);
    end
end

if isempty(range)
    iRowStart = 0;
    iColStart = 0;
    iRowEnd = size(A, 1)-1;
    iColEnd = size(A, 2)-1;
else
    iSeperator = strfind(range, ':');
    if isempty(iSeperator)
        cellStart = CellReference(range);
        iRowStart = cellStart.getRow();
        iColStart = cellStart.getCol();
        iRowEnd = iRowStart + size(A, 1)-1;
        iColEnd = iColStart + size(A, 2)-1;
    else
        cellStart = range(1:iSeperator-1);
        cellEnd = range(iSeperator+1:end);
        
        cellStart = CellReference(cellStart);
        cellEnd = CellReference(cellEnd);
        
        iRowStart = cellStart.getRow();
        iColStart = cellStart.getCol();
        iRowEnd = cellEnd.getRow();
        iColEnd = cellEnd.getCol();
    end
end

 
nRowA = size(A, 1)-1;
nColA = size(A, 2)-1;

 
if ~iscell(A)
    A = num2cell(A);
end

 
for iRow = iRowStart:iRowEnd
    currentRow = xlsSheet.getRow(iRow); 
    if isempty(currentRow)
        currentRow = xlsSheet.createRow(iRow);
    end
    
    for iCol = iColStart:iColEnd
        currentCell = currentRow.getCell(iCol);
        if isempty(currentCell)
            currentCell = currentRow.createCell(iCol);
        end
        
        if (iRow-iRowStart)<=nRowA && (iCol-iColStart)<=nColA
            data = A{iRow-iRowStart+1, iCol-iColStart+1};
            
            if ~isempty(data)          
                if isnumeric(data) && isnan(data)
                    data = '';
                end
                
                currentCell.setCellValue(data);
            end

        else
            currentCell.setCellErrorValue(FormulaError.NA.getCode());
        end
    end
end

fileOut = java.io.FileOutputStream(filename);
xlsWorkbook.write(fileOut);
fileOut.close();

status = 1;

end
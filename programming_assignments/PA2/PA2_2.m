[totalTax1, effectiveTaxRate1] = incomeTax(15000)
[totalTax2, effectiveTaxRate2] = incomeTax(44000)

function [totalTax, effectiveTaxRate] = incomeTax(Income)
%Write the commands for your function here.  

taxableIncome = Income - 6300 - 4000; %does the first deduction and exemption subtraction of the income

%use if statements and categorize each taxable income with its corresponding tax rate percentage
if taxableIncome >= 0.01 && taxableIncome <= 9225 %income range from 0.01 up to 9225
    totalTax = taxableIncome*0.1; %calculate the total tax from the taxable income
    effectiveTaxRate = totalTax / Income; % total tax divided by the income is the effective tax rate (as a decimal)
    
else if taxableIncome > 9225 && taxableIncome <= 37450%income range from 9225.01 up to 37450
        totalTax = (0.1 * 9225 + 0.15 * (taxableIncome - 9225));%calculate the total tax from the taxable income
        effectiveTaxRate = totalTax / Income;% total tax divided by the income is the effective tax rate (as a decimal)
        
else if taxableIncome > 37450 && taxableIncome <= 90750%income range from 37450.01 up to 90750
        totalTax = (0.1 * 9225 + 0.15 * (37450 - 9225) + 0.25 * (taxableIncome - 37450));%calculate the total tax from the taxable income
        effectiveTaxRate = totalTax / Income;% total tax divided by the income is the effective tax rate (as a decimal)
        
else if taxableIncome > 90750 && taxableIncome <= 189300;%income range from 90750.01 up to 189300
        totalTax = (0.1 * 9225 + 0.15 * (37450 - 9225) + 0.25 * (90750 - 37450) + 0.28 * (taxableIncome-90750));%calculate the total tax from the taxable income
        effectiveTaxRate = totalTax / Income;% total tax divided by the income is the effective tax rate (as a decimal)
        
else if taxableIncome > 189300 && taxableIncome <= 411500;%income range from 189300.01 up to 411500
        totalTax = (0.1 * 9225 + 0.15 * (37450 - 9225) + 0.25 * (90750 - 37450) + 0.28 * (189300-90750) + 0.33 * (taxableIncome - 189300));%calculate the total tax from the taxable income
        effectiveTaxRate = totalTax / Income;% total tax divided by the income is the effective tax rate (as a decimal)
        
else if taxableIncome > 411500 && taxableIncome <= 413200;%income range from 411500.01 up to 413200
        totalTax = (0.1 * 9225 + 0.15 * (37450 - 9225) + 0.25 * (90750 - 37450) + 0.28 * (189300-90750) + 0.33 * (411500 - 189300) + 0.35 * (taxableIncome - 411500));%calculate the total tax from the taxable income
        effectiveTaxRate = totalTax / Income;% total tax divided by the income is the effective tax rate (as a decimal)
        
        %income of 413200.01 and up
else totalTax = (0.1 * 9225 + 0.15 * (37450 - 9225) + 0.25 * (90750 - 37450) + 0.28 * (189300-90750) + 0.33 * (411500 - 189300) + 0.35 * (413200 - 411500) + 0.396 * (taxableIncome - 413200));%calculate the total tax from the taxable income
    effectiveTaxRate = totalTax / Income;% total tax divided by the income is the effective tax rate (as a decimal)
    
end  
end
end
end
end
end
end
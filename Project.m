close all; clc
choice = menu('Choices', 'Linear regression + models',...
    'Sol. for linear system','Sol. non-linear system',...
    'Numerical integration', 'Sol. ODE numericallly');

switch (choice)
    case 1 % Linear regression + models
        [x,y] = inputPoints();
        
        scatter(x,y);
        hold on;
        x_ = linspace(min(x), max(x));
        
        [a0,a1,r1] = linearRegression(x,y);
        y_linear = a0 + a1 * x_;
        plot(x_, y_linear);
        
        coeffl= ['The coefficient of correlation for linear regression is: ' num2str(r1) ];
        msgbox(coeffl);
        
        list = {'Exponential', 'Power', 'Growth-rate'};
        [model, tf] = listdlg('PromptString',{'Select a model.',...
            'You can select one model at a time!'},...
            'SelectionMode','single','ListString', list);
        
            switch (model)
            case 1 %Exponential
                [a,b,r2] = exponentialModel(x,y);
                y_exp = a * exp(b*x_);
                plot(x_, y_exp);
                coeffe= ['The coefficient of correlation of exponential model is: ' num2str(r2) ];
                msgbox(coeffe);
                legend('Original Points','Linear Model','Exponential Model');
            case 2 %Power
                [a,b,r3] = powerModel(x,y);
                y_power = a * (x_.^b);
                plot(x_, y_power);
                coeffp= ['The coefficient of correlation of power model is: ' num2str(r3) ];
                msgbox(coeffp);
                legend('Original Points','Linear Model','Power Model');
            case 3 %Growth-rate
                [a,b,r4] = growthRateModel(x,y);
                y_growth_rate = (a * x_)./(b + x_);
                plot(x_, y_growth_rate);
                coeffg= ['The coefficient of correlation of growth-rate model is: ' num2str(r4) ];
                msgbox(coeffg);
                legend('Original Points','Linear Model','Growth-rate Model');
                otherwise
            end
       
    case 2 % Sol. for linear system
        % Jacobi method
        
        rf=inputdlg({'Enter the number of rows of the coefficient matrix : ',...
            'Enter the number of columns of the coefficient matrix : '},...
            'inputs', [1 50], {'3','3'}, 'on');
        r=str2double(rf(1));
        c=str2double(rf(2));
        for k=1:r
            for m=1:c
                coeff= inputdlg(strcat('Enter the coefficient value in row ',...
                    num2str(k), 'column ', num2str(m), ': '), 'values');
                %3,-0.1,-0.2
                %0.1,7,-0.3
                %0.3,0.2,10
                A(k,m)=str2double(coeff);
            end
        end

        for l=1:r
             free= inputdlg(strcat('Enter the free term value in row ',...
                 num2str(l), ': '), 'value');
             %7.85
             %-19.3
             %71.4
             b(l)=str2double(free);
        end

        for p=1:c
             ini= inputdlg(strcat('Enter the initial value of variable no. ',...
                 num2str(p), ': '), 'value');
             %0
             %0
             %0
             initial(p)=str2double(ini);
        end
        
        lb = length(b);
        n_iterf = inputdlg('Enter the number of iterations : ',...
            'iterations', [1 50], {'3'}, 'on');
        n_iter= str2double(n_iterf);
        x = zeros(r, 1);
        for h=1:n_iter
            for i=1:c
                const = b(i)/A(i,i);
                for j=1 : r
                    if(i ~= j)
                    x(i) = A(i,j)*initial(j)/A(i,i);
                    x(i) = const - x(i);
                    end
                end
            end
        initial = x;
        end
        lx=length(x);
        for d=1:lx
            result =[strcat('The solution for variable ', num2str(d),...
                'is: ') num2str(x(d))];
            msgbox(result);
        end
    case 3 % Sol. non-linear system
        nonlinsolopts= menu('Methods', 'Newton-Raphson', 'Bisection');
        switch(nonlinsolopts)
            case 1 %Newton-Raphson
                % >>>>>>> example <<<<<<<<

                format long;
                ff= input('Enter a function of x : ', 's');
                % x.^2 + log(x) =0
                f = inline(ff);
%                 dft = diff(f);
%                 dft = functionalDerivative(f,x);
                dft=input('Enter the 1st derivative of the function x : ',...
                    's');
                % 2*x + 1/x
                df = inline(dft);
                
                inetialx_str = inputdlg('X0 = ', 'initial value', [1 50],...
                    {'0.5'}, 'on');
                initialx= str2double(inetialx_str);
                epsnon_str=inputdlg('maximum error bound = ','Max error',...
                    [1 50], {'0.0001'}, 'on');
                epsnon= str2double(epsnon_str);
                iterations = 1; 
                xnew = initialx - (f(initialx)/df(initialx));
                tolerance = abs(xnew - initialx); 
                while tolerance > epsnon %condition,,loop stops when epsnon > tolerance
                    disp = ([iterations,xnew])
                    iterations = iterations + 1 ;
                    initialx = xnew;
                    xnew = initialx - (f(initialx)/df(initialx));
                    tolerance = abs(xnew - initialx); 
                end
                    res1=['The iteration number is: ' num2str(iterations)];
                    msgbox(res1);
                    res2=['The solution is: ' num2str(xnew)];
                    msgbox(res2);
            case 2 %Bisection
                %takes the function from the user and return the root
                % >>>>>>> example <<<<<<<<
                % x^2 + log(x) =0
                infun  = input('Enter a function of x: ', 's');
                funct  = inline(infun);
                inarg  = inputdlg({'Enter the start of interval: ',...
                    'Enter the end of interval: '},...
                    'Initial arguments', [1 50], {'0.5','1'}, 'on');
                a      = str2double(inarg(1));
                b      = str2double(inarg(2));
                ques   = questdlg('would you enter a number of iterations or error?',...
                    'choose',... 
                    'Iterations', 'Error',[1 50]);
                switch ques
                    case 'Iterations'
                        itrs    = inputdlg('Enter the number of iterations: ',...
                        'Number of itrations', [1 50],{'4'}, 'on');
                        itr     = str2double(itrs);
                    case 'Error'
                        errors  = inputdlg('Enter the maximum error bound: ',...
                        'Error bound', [1 50],{'0.03125'}, 'on');
                        er      = str2double(errors);
                        itr     = abs(log2((b-a)/er));
                    otherwise
                end

                if funct(a)*funct(b)<0
                else
                    errormsg = questdlg('Please enter valid guesses!\n');
                    inarg= inputdlg({'Enter the start of interval: ',...
                    'Enter the end of interval: '},...
                    'Initial arguments', [1 50],{'0.5','1'},'on');
                    a= str2double(inarg(1));
                    b= str2double(inarg(2));
                end
                for i=0:itr
                        c= (a+b)/2;
                        if funct(a)*funct(c)<0
                            b=c;
                        else
                            a=c;
                        end
                        if funct(b)*funct(c)<0
                            a=c;
                        else
                            b=c;
                        end
                end
                result=['The required root is: ' num2str(c)];
                msgbox(result);
            otherwise
        end
      
    case 4 %Numerical integration
        intopts= menu('Methods', 'Trapezoidal', 'Simpson’s 1/3 rule');
        switch (intopts)
            case 1 % Trapezoidal
                %the user will enter function y(x)=function of x
                %the user will enter the integration limits 
                %the user will enter number of segments for higher accuracy
                %then we will apply the trapezoidal rule 
                % >>>>>>>> example <<<<<<<<<
                % ingration 

                %interval [0,1] 
                %5 segments
                yf = input('enter a function of x:  ', 's');
                                % x*exp(-x) 
                y= inline(yf);
                inputs= inputdlg({'enter the start ','enter the end ',...
                    'Enter number of segments '}, 'Input values',...
                    [1 50], {'0','1','5'}, 'on');
                starttime = str2double(inputs(1));
                endtime = str2double(inputs(2));
                segments= str2double(inputs(3));
                %calculate the step size 
                step = (endtime-starttime)/segments;

                sum1 = 0;
                for i=1 : segments-1 
                    f = starttime + i*step;

                    sum1 = sum1 + y(f);
                end
                result=(step/2)*(y(starttime)+y(endtime)+2*sum1);
                resulto=['The result is: ' num2str(result)];
                msgbox(resulto);
            case 2 % Simpson’s 1/3 rule
                % >>>>>>> example <<<<<<
                % x*cos(exp(x)) [0,2] h=0.25
                yf = input('enter a function of x :  ', 's');
                y=inline(yf);
                inputs= inputdlg({'enter the start: ','enter the end: ',...
                    'Enter step size (h): '}, 'Input values', [1 50],...
                    {'0','2','0.25'}, 'on');
                starttime = str2double(inputs(1));
                endtime = str2double(inputs(2));
                segments= (endtime-starttime)/str2double(inputs(3));
                %calculate the step size 
                step = str2double(inputs(3));
                s_odd=0;  s_even =0;

                %to get sum for the odd part
                for i = 1 : 2 :segments-1
                    s_odd = s_odd + 4* y(starttime + i*step);
                end

                %to get sum for the even part
                for i = 2 : 2 :segments-2
                    s_even = s_even+ 2* y(starttime + i*step);
                end
                result=(step/3)*(y(starttime)+y(endtime)+s_even+s_odd);
                res=['The result is: ' num2str(result)];
                msgbox(res);
            otherwise
        end
    case 5 %Sol. ODE numericallly
        odeopts = menu('Methods', 'Euler’s' , 'Heun’s');
        switch (odeopts)
            case 1 % Euler’s
                % >>>>>>> example <<<<<<
                % dy/dx = y*(x^3)-1.5*y  y(0)=1   x=1  h=0.25
                fun=input('Enter the function of x and y: ','s');
                %y*(x^3)-1.5*y
                f = inline(fun,'x','y');

                inps= inputdlg({'Enter the step size: ',...
                    'Enter the start point of x: ',...
                    'Enter the end point of x: ',...
                    'Enter initial value of y: ' },'Inputs', [1 50],...
                    {'0.25','0','1','1'},'on');

                m=str2double(inps(1));
                a=str2double(inps(2));
                b = str2double(inps(3));
                ya =str2double(inps(4));

                h=(b-a)/m;
                y=zeros(h+1,1);
                x=(a:m:b);
                y(1)=ya;

                for j=1:h
                    y(j+1)=y(j)+m*f(x(j),y(j));
                end

                plot(x,y);
            case 2 % Heun’s
                 % >>>>>>> example <<<<<<
                % dy/dx = y*(x^3)-1.5*y  y(0)=1   x=1  h=0.25
                
                fun=input('Enter the function of x and y: ','s');
                %y*(x^3)-1.5*y
                f = inline(fun,'x','y');

                inps= inputdlg({'Enter the step size: ',...
                    'Enter the start point of x: ',...
                   'Enter the end point of x: ',...
                   'Enter initial value of y: ' },'Inputs', [1 50],...
                   {'0.25','0','1','1'},'on');

                m=str2double(inps(1));
                a=str2double(inps(2));
                b = str2double(inps(3));
                ya =str2double(inps(4));

                h=(b-a)/m;  %NO. of steps
                y=zeros(h+1,1);
                x=(a:m:b);
                y(1)=ya;
                y0=zeros(h+1,1);


                for j=1:h

                    y0(j+1)=y(j)+m*f(x(j),y(j));  %predictor
                    y(j+1)=y0(j);
                    y0(1)=ya;
                    y(j+1)=y(j)+m*(f(x(j),y(j))+f(x(j+1),y0(j+1)))/2;  %corrector
                end
                figure
                plot(x,y);
            otherwise
        end
    otherwise
end

%%%%%%%%%%%%%%%%%%%functions definition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = inputPoints()
    nstr = inputdlg('Enter the number of points: ', 'No. of points', [1 50]);
    n=str2double(nstr);
    x= zeros(1,n);
    y = zeros(1,n);

    i=1;

    while i<=numel(x)
        inptx = inputdlg(strcat('Enter the point for the x number(',num2str(i),'): '), 'x');
        x(i) = str2double(inptx);
        i=i+1;
    end
    
    i=1;
    
    while i<=numel(y)
        inpty = inputdlg(strcat('Enter the corresponding point of y for the x number(',num2str(i),'): '), 'y');
        y(i) = str2double(inpty);
        i=i+1;
    end
end

function [a0,a1,r] = linearRegression(x,y)
    n = numel(x);
    sigmaX = sum(x, "all");
    sigmaY = sum(y,"all");
    sigmaX_squared = sum(x.^2,"all");
    sigmaXY = sum(x.*y, "all");
    
    A=[n sigmaX; sigmaX, sigmaX_squared];
    B=[sigmaY; sigmaXY];
    
    X=linsolve(A, B);
    a0=X(1);
    a1=X(2);
    r=correlationCoefficient(x,y,a0,a1);
end

function [a,b,r] = exponentialModel(x,y)
    Y = log(y);
    [a0,a1] = linearRegression(x,Y);
    a=exp(a0);
    b=a1;
    r=correlationCoefficient(x,Y,a0,a1);
end

function [a,b,r] = powerModel(x,y)
    Y = log(y);
    X = log(x);
    [a0,a1] = linearRegression(X,Y);
    a=exp(a0);
    b=a1;
    r=correlationCoefficient(X,Y,a0,a1);
end

function [a,b,r] = growthRateModel(x,y)
    Y = 1./y;
    X = 1./x;
    [a0,a1] = linearRegression(X,Y);
    a=1/a0;
    b=a*a1;
    r=correlationCoefficient(X,Y,a0,a1);
end

function r = correlationCoefficient(x,y,a0,a1)
    y_ = mean(y);

    St = 0;
    Sr = 0;
    for n = 1:length(y)
        St = St + (y(n) - y_)^2;
        Sr = Sr + (y(n) - a0 - a1 * x(n))^2;
    end

    r = sqrt((St-Sr)/St);
end


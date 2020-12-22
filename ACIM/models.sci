// ------- Public Functions --------
exec("fn_combine.sci");

function y = evalModel(model, evalData)
    y = model.Function(evalData, model.Params);
    // exclude negative numbers (e.g. there can´t be a negative number 
    // for number of leaves etc.)
    if y < 0 then
        y = 0;
    end
endfunction

function model = trainModel(model, data)
    // adjust data and functions (transpose) for 
    global sourceData
    sourceData = data(:,1:$-1);
    deff('[e] = errFct(p,data)', ...
         [
          // transpose
         'data = data''',
          // transpose
         'p = p''',
         'e = model.ErrorFunction(data,p)',
         // rate all negative model data with a high error value
         //   for the training data
         'e(model.Function(data(:,1:$-1),p)<0) = 1e20',
         //   for interpolated data
         'global sourceData'
         'iData = linspace(min(sourceData,''r'')'',max(sourceData,''r'')'',25)''',
         'if (~isempty(find(model.Function(iData,p)<0))) e(1) = 1e20 end'
         ]);
    initParams = model.InitialParams';
    dat = data';
    
    disp("start  for model " + model.Name)
    if (isfield(model, "Binf")) then
        [p,err] = datafit(errFct, dat, 'b', model.Binf, model.Bsup, initParams);
    else
        [p,err] = datafit(errFct, dat, initParams);
    end
    disp("  ready")

    // Save fitting
    model.Params = p';
    
    model.MeanSquareError = err / size(data, 1);
    
    numDataP = size(data, 1); // number of data points
    model.AIC = calcAIC(model.MeanSquareError, length(model.Params), numDataP);
    
    targetMean = mean(data(:,$));
    // sum of squared errors of baseline model
    SST = sum((data(:,$)-targetMean).^2);
    // Coefficient of determination
    model.R2 = 1 - (err/SST);

    // Coefficient of determination adjusted for multiple regression
    // Degrees of freedom considered, definition from: http://www.tf.uns.ac.rs/~omorr/radovan_omorjan_003_prII/s_examples/Scilab/Gilberto/scilab17.pdf
    numVar = model.InputCount; // number of influencing variables
    model.R2adj = 1 - ((err/(numDataP-numVar-1))/(SST/(numDataP-1)));
endfunction

function allModels = getAllModels(baseModels, combineCount)
    allModels = list();
    for i = 1:combineCount
        modelCombinations = getAllCombinations(baseModels, i);
        
//        disp("combine level " + string(i) + " has " + string(length(modelCombinations)) + " models")
//        showModelCombinations(modelCombinations)
        
        models = list();
        for modelCombination = modelCombinations
            models($+1) = createCombinedModel(modelCombination);
        end
        allModels($+1) = models;
    end
endfunction

function baseModels = getBaseModels()
    baseModels = list();
    baseModels($+1) = LinModel();
    baseModels($+1) = QuadModel();
    baseModels($+1) = SigmModel();
    //baseModels($+1) = RatModel();
    baseModels($+1) = GaussModel();
    baseModels($+1) = ExpDcModel();
    baseModels($+1) = ExpGrModel();
    baseModels($+1) = HyperbolaModel();
endfunction

function aic = calcAIC(meanSquareError, paramCount, observationCount)
    aic = 2 * paramCount + observationCount * log(meanSquareError);
end

// ------- Private Functions --------

// ----- Debug Functions ------

function showModelCombinations(combinedModels)
    for i = 1:length(combinedModels)
        disp("combined model " + string(i))
        disp("  model count: " + string(length(combinedModels(i))))
        for model = combinedModels(i)
            disp("  " + model.Name)
        end
    end
endfunction

function showModel(model)
    disp("model name: " + model.Name)
    disp("   count: " + string(model.InputCount))
    disp("  param count: " + string(length(model.InitialParams)))
    disp(["  initial params: " string(model.InitialParams)])
    if (isfield(model, "Params")) then
        disp(["  params: " string(model.Params)])
        disp("  mean square error: " + string(model.MeanSquareError))
        disp("  aic: " + string(model.AIC))
        disp("  sample : " + string(size(model.Sources(1).Data, 1)))
    else
        dummyData = createDummyData(model.InputCount);
        sampleError = model.ErrorFunction(dummyData, model.InitialParams);
        disp("  dummy error: " + string(sampleError))
    end
endfunction


function showModels(models)
    for model = models
        if (model <> [])
            showModel(model);
        end
    end
endfunction

function dataSample = createDummyData(Count)
    dataSample = ones(1,Count+1);
endfunction

// ----- Model Combination Functions ------

function model = createCombinedModel(models)
    modelCount = length(models);

    select modelCount,
        case 0 then
            model = [],
        case 1 then 
            model = models(1)
        else            
            // Name
            model.Name = models(1).Name;
            for i = 2:modelCount
                model.Name = model.Name + "_" + models(i).Name;
            end
            
            //  Count
            model.InputCount = 0;            
            for i = 1:modelCount
                model.InputCount = model.InputCount + models(i).InputCount;
            end
            
            // Params
            paramCount = modelCount + 1;
            for i = 1:modelCount
                paramCount = paramCount + length(models(i).InitialParams);
            end
            model.InitialParams = createInitialParams(paramCount);
            
            // param constrains 
//            disp("param contrains")
            model.Binf = [];
            model.Bsup = [];
            for i = 1:modelCount
                model.Binf = [model.Binf, models(i).Binf];
                model.Bsup = [model.Bsup, models(i).Bsup];
            end
            
            // for weights
            model.Binf = [model.Binf, zeros(1,modelCount+1)];//lower bound
            model.Bsup = [model.Bsup, ones(1,modelCount+1)]; // upper bound
//            disp(model.Binf)
//            disp(model.Bsup)
//            disp("Binf ready")    
            // Function
            model.Function = createCombinedFunction(model.Name, models, 'y');
            
            // Error Function
            model.ErrorFunction = createErrorFunction(model.Name);
            
            showModel(model);
    end
endfunction

function fct = createCombinedFunction(modelName, models, output)
    functionBody = [];
    // Make model functions available
    for model = models
        functionBody = [functionBody "global " + getFunctionName(model.Name)];
    end
    // Concat addends of each model function
    // e.g.: y = anyXFct(x(1:1,:),p(1:3)).*p(4) + anyYFct(x(2:2,:),p(5:6)).*p(7);
    //   Collect model terms (addends)
    xIndex = 1;
    pIndex = 1;
    modelTerms = list();
    for model = models
        xFrom = string(xIndex);
        xTo = string(xIndex+model.InputCount-1);
        pFrom = string(pIndex);
        pTo = string(pIndex+length(model.InitialParams)-1);
        addend = [getFunctionName(model.Name) "( x(:," xFrom ":" xTo "), p(" pFrom ":" pTo ") )"];
        modelTerms($+1) = addend;
        
        xIndex = xIndex + model.InputCount;
        pIndex = pIndex + length(model.InitialParams);
    end
    
    //   Use terms to create function - first linear independent, then all linear dependent
    functionString = ["y = "];
    for modelTerm = modelTerms
        functionString = [functionString modelTerm " .* p(" string(pIndex) ")" " + "];
        pIndex = pIndex + 1;
    end
    for modelTerm = modelTerms
        functionString = [functionString modelTerm " .* "];
    end
    functionString = [functionString "p(" string(pIndex) ");"];
    //disp(strcat(functionString))
    functionBody = [functionBody strcat(functionString)];
    //disp(functionBody)
    
    fct = createFunction(modelName, output, functionBody);
endfunction

function combinations = getCombinations(values, count)
    combinations = list();
    
    valuesVec = list2vec(values);
    combinationsIndices = strtod(combine((1:length(values))',count,%F))(:,1:$-1);
    for combinationIndices = combinationsIndices'
        combinations($+1) = (valuesVec(combinationIndices))';
    end
endfunction

function combinations = getAllCombinations(values, depth)
    combinations = list();
    if (depth > 0) then
        lowerCombinations = getAllCombinations(values, depth - 1);
        if (length(lowerCombinations) == 0)
            for value = values
                combinations($+1) = list(value);
            end
        else
            for value = values
                for lowerCombination = lowerCombinations
                    combinations($+1) = lstcat(list(value),lowerCombination);
                end
            end
        end
    end
endfunction

// ----- Model Creation Functions ------

function name = getFunctionName(modelName)
    name = modelName + "_function";
endfunction

function name = getErrorFunctionName(modelName)
    name = modelName + "_error_function";
endfunction

function fct = createFunction(modelName, output, functionBody)
    functionName = getFunctionName(modelName);
    execstr('global ' + functionName);
    
    deff('[' + output + '] = ' + functionName + '(x,p)', functionBody);
    execstr('fct = ' + functionName);
endfunction

function errFct = createErrorFunction(modelName)
    functionName = getFunctionName(modelName);
    errorFunctionName =getErrorFunctionName(modelName);
    execstr('global ' + errorFunctionName);
    
    deff('[e]=' + errorFunctionName + '(data,p)', ...
         [
         'global ' + functionName,
         'mData = ' + functionName + '(data(:,1:$-1),p)',
         'e = data(:,$)-mData',
         ]);
    execstr('errFct = ' + errorFunctionName);
endfunction

function params = createInitialParams(paramCount)
    params = ones(1, paramCount);
endfunction

// ------- Base Model Functions --------

function m=LinModel()
    m.InputCount = 1;
    m.Name = "linear";
    m.InitialParams = createInitialParams(3);
    m.Binf = [-%inf,-%inf,0];
    m.Bsup = [%inf,%inf,%inf];
    m.Function = createFunction(m.Name, 'y', 'y =(p(1).*(x+p(2))+p(3))');
    m.ErrorFunction = createErrorFunction(m.Name);
endfunction

function m=QuadModel()
    m.InputCount = 1;
    m.Name = "quadratic";
    m.InitialParams = createInitialParams(4); 
    m.Binf = [-%inf,-%inf,-%inf,-%inf];
    m.Bsup = [%inf,%inf,%inf,%inf];
    m.Function = createFunction(m.Name, 'y', 'y = (p(1).*(x+p(2)).^2+p(3).*(x+p(2))+p(4))');
    m.ErrorFunction = createErrorFunction(m.Name);
endfunction

function m=SigmModel()
    m.InputCount = 1;
    m.Name = "sigmoidal";
    m.InitialParams = createInitialParams(4);
    m.Binf = [0,-%inf,-%inf,-%inf];
    m.Bsup = [1,%inf,%inf,%inf];
    m.Function = createFunction(m.Name, 'y', 'y=((p(1).*exp(p(2).*(x+p(3))))./(p(1).*exp(p(2).*(x+p(3)))+(1-p(1))))+p(4)');
    m.ErrorFunction = createErrorFunction(m.Name);
endfunction

function m=HyperbolaModel()
    m.InputCount = 1;
    m.Name = "hyperbola";
    m.InitialParams = createInitialParams(3);
    m.Binf = [-%inf,-%inf,-%inf];
    m.Bsup = [%inf,%inf,%inf]; 
    m.Function = createFunction(m.Name, 'y', 'y = (p(1)./(x+p(2)))+p(3)');
    m.ErrorFunction = createErrorFunction(m.Name);
endfunction
//
//function m=RatModel()
//    m.Count = 1;
//    m.Name = "rational";
//    m.InitialParams = createInitialParams(4);
//    m.Function = createFunction(m.Name, 'y', 'y=(p(1).*x+p(3))./(p(2).*x+p(4))');
//    m.ErrorFunction = createErrorFunction(m.Name);
//endfunction
//
function m=GaussModel()
    m.InputCount = 1;
    m.Name = "gaussian";
    m.InitialParams = createInitialParams(4);
    m.Binf = [-%inf,-%inf,-%inf,-%inf]; // check here later if constraints are ok
    m.Bsup = [%inf,%inf,%inf,%inf]; // check here later if constraints are ok
     //function from Wikipedia
    m.Function = createFunction(m.Name, 'y', 'y = (p(4)+(p(1)*exp(-((x-p(2)).^2/(2*p(3))).^2)))');
    m.ErrorFunction = createErrorFunction(m.Name);
endfunction

function m=ExpDcModel()
    m.InputCount = 1;
    m.Name = "exponentialDC"; // p3 min 0 max x 
    m.InitialParams = [1,-1,1,1];//createInitialParams(3);
    m.Binf = [-%inf,-%inf,0,-%inf];
    m.Bsup = [%inf,0,%inf,%inf];
   // 2019 DELETE LATER: m.Function = createFunction(m.Name, 'y', 'y= p(3)+(p(1)*exp(p(2)*x))');
    m.Function = createFunction(m.Name, 'y', 'y= p(3)+(p(1)*exp(p(2)*(x+p(4))))');
    // alt 'y=((p(1).*exp(p(2).*(x+p(3))))+p(4))');
    m.ErrorFunction = createErrorFunction(m.Name);
endfunction

function m=ExpGrModel()
    m.InputCount = 1;
    m.Name = "exponentialGr"; // p3 min 0 max x 
    m.InitialParams = [1,1,1,1];//createInitialParams(3);
    m.Binf = [-%inf,0,-%inf,-%inf];
    m.Bsup = [%inf,%inf,%inf,%inf];
    // 2019 DELETE LATER: m.Function = createFunction(m.Name, 'y', 'y= p(3)+(p(1)*exp(p(2)*x))');
    m.Function = createFunction(m.Name, 'y', 'y= p(3)+(p(1)*exp(p(2)*(x+p(4))))'); 
    
    // alt 'y=((p(1).*exp(p(2).*(x+p(3))))+p(4))');
    m.ErrorFunction = createErrorFunction(m.Name);
endfunction

// ------- Test Methods --------

function testLinearModel()
    data = [1,20;
            2,31;
            3,39;
            4,45];
    
    untrainedModel = LinModel();
    showModel(untrainedModel);
    trainedModel = trainModel(untrainedModel, data);
    disp(trainedModel)
    showModel(trainedModel);
endfunction

function testModel()
    path = '.\Testrelationship\';
    //untrainedModel = ExpGrowModel()
    Rels = createAllModels(path);
    //data = list()
    baseModels = getBaseModels()

    ModelType = "Test";
    for i=1:length(Rels)
        for k=1:length(baseModels)
            trainedModel = list();
            untrainedModel = baseModels(k);
            data = [Rels(i).Models(1).Sources.Data,Rels(i).Models(1).Target.Data]
            trainedModel = trainModel(untrainedModel, data);
            trainedModel.Sources.Data = Rels(i).Models(1).Sources.Data;
            trainedModel.Sources.Name = Rels(i).Models(1).Sources.Name;
            trainedModel.Sources.Unit = Rels(i).Models(1).Sources.Unit;
            trainedModel.Target.Data = Rels(i).Models(1).Target.Data;
            trainedModel.Target.Name = Rels(i).Models(1).Target.Name;
            trainedModel.Target.Unit = Rels(i).Models(1).Target.Unit;
            trainedModel.Publication = " ";
            figure(k)
            displayGraphForModel(trainedModel,ModelType)
        end
//        for j=1:length(Rels(i).Models)
//            disp(Rels(i).Models(j).Name)
//            displayGraphForModel(Rels(i).Models(j),ModelType)
//        end
    end
    //showModel(trainedModel);
endfunction

function testExpDcModel()
    path = '.\Testrelationship\';
    //untrainedModel = ExpGrowModel()
    Rels = createAllModels(path);
    //data = list()
    
    ModelType = "Test";
    for i=1:length(Rels)
        trainedModel = list();
        untrainedModel = ExpDcModel();
        data = [Rels(i).Models(1).Sources.Data,Rels(i).Models(1).Target.Data]
        trainedModel = trainModel(untrainedModel, data);
        trainedModel.Sources.Data = Rels(i).Models(1).Sources.Data;
        trainedModel.Sources.Name = Rels(i).Models(1).Sources.Name;
        trainedModel.Sources.Unit = Rels(i).Models(1).Sources.Unit;
        trainedModel.Target.Data = Rels(i).Models(1).Target.Data;
        trainedModel.Target.Name = Rels(i).Models(1).Target.Name;
        trainedModel.Target.Unit = Rels(i).Models(1).Target.Unit;
        trainedModel.Publication = " ";
        figure(i)
        displayGraphForModel(trainedModel,ModelType)
    end
    //showModel(trainedModel);
endfunction

function testExpGrModel()
    path = '.\Testrelationship\';
    //untrainedModel = ExpGrowModel()
    Rels = createAllModels(path);
    //data = list()
    
    ModelType = "Test";
    for i=1:length(Rels)
        trainedModel = list();
        untrainedModel = ExpGrModel();
        data = [Rels(i).Models(1).Sources.Data,Rels(i).Models(1).Target.Data]
        trainedModel = trainModel(untrainedModel, data);
        trainedModel.Sources.Data = Rels(i).Models(1).Sources.Data;
        trainedModel.Sources.Name = Rels(i).Models(1).Sources.Name;
        trainedModel.Sources.Unit = Rels(i).Models(1).Sources.Unit;
        trainedModel.Target.Data = Rels(i).Models(1).Target.Data;
        trainedModel.Target.Name = Rels(i).Models(1).Target.Name;
        trainedModel.Target.Unit = Rels(i).Models(1).Target.Unit;
        trainedModel.Publication = " ";
        figure(i)
        displayGraphForModel(trainedModel,ModelType)
    end
    //showModel(trainedModel);
endfunction

function testLinearLinearModel()
    data = [1,2,20;
            2,2,31;
            3,2,39;
            4,2,45];
    
    untrainedModel = createCombinedModel([LinModel() LinModel()]);
    showModel(untrainedModel);
    trainedModel = trainModel(untrainedModel, data);
    showModel(trainedModel);
endfunction

function testAllModels()
    combineCount = 3;
    
    allModels = getAllModels(getBaseModels(), combineCount);
    
    for i = 1:length(allModels)
        disp("model combination count: " + string(i))
        showModels(allModels(i))
    end
endfunction

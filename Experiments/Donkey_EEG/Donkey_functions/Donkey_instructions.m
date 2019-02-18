%%%%%%%%%%%%%%%%
%% Instructions
%%%%%%%%%%%%%%%%

% Set font to standard
Screen('TextFont', w, standardFont);

switch instructions
        
    %% Main instructions
    case 'introduction'
        
        % Welcome screen
        Screen(w,'FillRect',col.background);   % background colour
        Screen('TextStyle',w, 1); % bold
        Screen('TextSize', w , titleSize);
        
        if whichPhase == 1
            Screen(w,'DrawText',['Part ' num2str(whichPhase) ': A visit to the farm.'],leftMargin,screenYpixels * 0.1,col.white);            
        elseif whichPhase == 2
            Screen(w,'DrawText',['Part ' num2str(whichPhase) ': Working on the farm.'],leftMargin,screenYpixels * 0.1,col.white);
        elseif whichPhase == 3
            Screen(w,'DrawText',['Part ' num2str(whichPhase) ': Your escape back home.'],leftMargin,screenYpixels * 0.1,col.white);
        end
        
        Screen('TextSize', w,textSize);
        Screen('TextStyle',w, 0); %normal
        Screen(w,'DrawText','Press spacebar to continue ...',leftMargin,screenYpixels * 0.9,col.white);
        
        Screen(w,'Flip');        
        KbWait;
        while KbCheck; end;
        
        % //Introduction page 1
        mytxtfile = ['Introduction_' num2str(whichPhase) '_p1.m'];
        mytext = printText(mytxtfile);     
        DrawFormattedText(w, mytext, leftMargin, topMargin, col.white,[],[],[],1.5);
               
        if whichPhase == 1
            Screen(w,'DrawText','Press spacebar to continue ...',leftMargin,screenYpixels * 0.9,col.white);
            Screen(w,'Flip');
            KbWait;
            while KbCheck; end;
            
            % //Introduction page 2
            mytxtfile = ['Introduction_' num2str(whichPhase) '_p2.m'];
            mytext = printText(mytxtfile);
            DrawFormattedText(w, mytext, leftMargin, topMargin, col.white,[],[],[],1.5);
        end
        
        if whichPhase == 2          
            Screen(w,'DrawText','Press spacebar to continue ...',leftMargin,screenYpixels * 0.9,col.white);           
            Screen(w,'Flip');
            KbWait;
            while KbCheck; end;
            
            % //Conclusion
            mytxtfile = ['Conclusion_' num2str(whichPhase) '.m'];
            mytext = printText(mytxtfile); % adjust last input to number of lines in the text file
            DrawFormattedText(w, mytext, leftMargin, topMargin, col.white,[],[],[],1.5);
        end
        
        %% Instructions at start of each block
    case 'startblock'
        
        Screen(w,'FillRect',col.background);
        Screen('TextSize', w,titleSize);
        Screen('TextStyle',w,1); % bold
        
        if whichPhase == 1
            line1 = ['Ready?'];
        elseif whichPhase == 2
            line1 = ['Hour ' num2str(data.block(t))];
        elseif whichPhase == 3
            line1 = ['Round ' num2str(data.block(t))];
        end
        DrawFormattedText(w, line1,leftMargin, screenYpixels * 0.1, col.white,[],[],[],1.5);
        
        %% Instructions at the end of each block
    case 'endblock'
        
        Screen(w,'FillRect',col.background);
        Screen('TextSize', w,titleSize);
        Screen('TextStyle',w,1); % bold
        line1 = 'Round completed';
        DrawFormattedText(w, line1,leftMargin, screenYpixels * 0.1, col.white,[],[],[],1.5);
        
        Screen('TextSize', w,textSize);
        Screen('TextStyle',w, 0); % normal
        
        line1 = 'Take a short break.';
        if whichPhase == 3
            line2 = ['\nCoins earned in this block: ' num2str(data.blockReward(t))];
            DrawFormattedText(w,[line1 line2], 'center', 'center', col.white,[],[],[],1.5);
        else
            DrawFormattedText(w,[line1], 'center', 'center', col.white,[],[],[],1.5);
        end
        
        %% Instructions at end of experiment
    case 'endexp'
        
        Screen(w,'FillRect',col.background);
        Screen('TextSize', w,titleSize);
        Screen('TextStyle',w,1); % bold
        line1 = ['Part ' num2str(whichPhase) ' finished!'];
        DrawFormattedText(w, line1,'center', screenYpixels * 0.25, col.white,[],[],[],1.5);
        
        if whichPhase == 3
            Screen('TextSize', w,textSize);
            Screen('TextStyle',w, 0); % normal
            line1 = ['\nYou''ve earned ' num2str(money) ' pounds extra.'];
            line2 = ['\nThat''s enough money to go back home!'];
            DrawFormattedText(w,[line1 line2], 'center', 'center', col.white,[],[],[],1.5);     
        else
            Screen('TextSize', w,textSize);
            Screen('TextStyle',w, 0); % normal
            line1 = ['\nNext, the screen will shut down and open the next part.'];
            DrawFormattedText(w,[line1], 'center', 'center', col.white,[],[],[],1.5);
        end
end

Screen('TextSize', w,textSize);
Screen('TextStyle',w, 0); % normal
Screen(w,'DrawText','Press spacebar to continue ...',leftMargin,screenYpixels * 0.9,col.white);
Screen(w,'Flip');

KbWait;
while KbCheck; end;

Screen(w,'FillRect',col.background);
Screen(w,'Flip');
WaitSecs(.5);
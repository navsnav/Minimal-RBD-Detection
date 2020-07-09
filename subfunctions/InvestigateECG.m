rem_idx = ismember(Sleep_table.AnnotatedSleepStage,5);
rbd_idx = ismember(Sleep_table.SubjectCondition,5);

figure;
boxplot(Sleep_table.EMG_AtoniaIndex(rem_idx),Sleep_table.SubjectCondition(rem_idx));

RBD_table = Calculate_EMG_Values_table(Sleep_table);

%%
subject_names = fieldnames(Sleep_Struct);
id = 38;
sub_idx = ismember(Sleep_table.SubjectIndex,id);

figure;
a(1) = subplot(2,1,1);
stairs(Sleep_table.AnnotatedSleepStage(sub_idx))
a(2) = subplot(2,1,2);
stairs(Sleep_table.SDNN(sub_idx))
hold on;
stairs(Sleep_table.SDNN_150s(sub_idx))
linkaxes(a,'x');


figure;
a(1) = subplot(2,1,1);
stairs(Sleep_table.AnnotatedSleepStage(sub_idx))
a(2) = subplot(2,1,2);
stairs(Sleep_table.LFHF(sub_idx))
hold on;
stairs(Sleep_table.LFHF_150s(sub_idx))
linkaxes(a,'x');


sub_idx = ismember(Sleep_table.SubjectCondition,0);


figure;
boxplot(Sleep_table.RR(sub_idx),Sleep_table.AnnotatedSleepStage(sub_idx));

figure;
boxplot(Sleep_table.RR_150s(sub_idx),Sleep_table.AnnotatedSleepStage(sub_idx));


sub_idx = ismember(Sleep_table.SubjectCondition,5);

figure;
boxplot(Sleep_table.SampEn(sub_idx),Sleep_table.AnnotatedSleepStage(sub_idx));

figure;
boxplot(Sleep_table.SampEn_150s(sub_idx),Sleep_table.AnnotatedSleepStage(sub_idx));

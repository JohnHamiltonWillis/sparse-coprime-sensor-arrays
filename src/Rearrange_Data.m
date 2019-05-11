function [rearranged_data, correct_order] = Rearrange_Data(totalData)
% Rearrange_Data Rearrange data from LaTech Senior Design array so that the
% data is linear from left to right
% second output argument is the correct_order vector. 
correct_order = [16	32	15	31	14	30	13	29	12	28	11	27	10	26	9	25	8	24	7	23	6	22	5	21	4	20	3	19	2	18	1	17	33	49	34	50	35	51	36	52	37	53	38	54	39	55	40	56	41	57	42	58	43	59	44	60	45	61	46	62	47	63	48	64];

rearranged_data = totalData(:,correct_order);
end

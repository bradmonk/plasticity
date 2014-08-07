function [fntname] = goodfont(fnum)

switch fnum
    case 1
		fntname = 'Euclid Extra';
	case 2
		fntname = 'Arial Unicode MS';
	case 3
		fntname = 'Times New Roman';
	case 4
		fntname = 'Helvetica';
	case 5
		fntname = 'Century Gothic';
    case 6
		fntname = 'Courier';
	case 7
		fntname = 'Arial Black';
	case 8
		fntname = 'Glowworm';
	case 9
		fntname = 'Dialog';
    case 10
		fntname = 'DIN Alternate';
	case 11
		fntname = 'Diwan Kufi';
	case 12
		fntname = 'Baghdad';
	case 13
		fntname = 'Euclid Math One';
	case 14
		fntname = 'Geometr885 BT';
	case 15
		fntname = 'BlairMdITC TT';
	case 16
		fntname = 'Hei';
	case 17
		fntname = 'Silom';
	case 18
		fntname = 'Avenir';
	case 19
		fntname = 'Yuanti SC';
otherwise
        warning('Only 19 good fonts!');
end




end
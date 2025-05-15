function toggle(talk) {
    var c = talk.className;
    if (c.match(" hidden")) {
	talk.className = c.replace(" hidden", "");
    } else {
	talk.className = c + " hidden";
    }
}

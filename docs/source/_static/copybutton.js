document.addEventListener("DOMContentLoaded", () => {
  const codeBlocks = document.querySelectorAll("div.highlight");

  codeBlocks.forEach((block) => {
    if (block.querySelector(".copy-code-button")) {
      return;
    }

    const button = document.createElement("button");
    button.type = "button";
    button.className = "copy-code-button";
    button.setAttribute("aria-label", "Copy code block");
    button.textContent = "Copy";

    button.addEventListener("click", async () => {
      const code = block.querySelector("pre");
      if (!code) {
        return;
      }

      const lines = code.innerText.split("\n");
      const isDoctestBlock = lines.some(
        (line) => line.startsWith(">>> ") || line.startsWith("... ")
      );

      const text = isDoctestBlock
        ? lines
            .filter(
              (line) => line.startsWith(">>> ") || line.startsWith("... ")
            )
            .map((line) => line.slice(4))
            .join("\n")
        : lines.join("\n");

      try {
        await navigator.clipboard.writeText(text);
        button.textContent = "Copied";
        window.setTimeout(() => {
          button.textContent = "Copy";
        }, 1500);
      } catch (err) {
        button.textContent = "Failed";
        window.setTimeout(() => {
          button.textContent = "Copy";
        }, 1500);
      }
    });

    block.classList.add("has-copy-button");
    block.appendChild(button);
  });
});

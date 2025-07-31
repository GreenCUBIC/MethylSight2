import type { RequestHandler } from "./$types";
import { text } from "@sveltejs/kit";

export const POST: RequestHandler = async ({ request }) => {
  const formData = await request.formData();

  const sequenceText = formData.get("sequenceText") as string;
  const file = formData.get("file") as File | null;

  let sequence = "";

  // Use text area input if provided
  if (sequenceText && sequenceText.trim() !== "") {
    sequence = sequenceText.trim().replace("\r\n", "\n");
  }
  // Otherwise, try to read the file
  else if (file && file.size > 0) {
    const buffer = await file.arrayBuffer();
    sequence = new TextDecoder().decode(buffer).trim().replace("\r\n", "\n");
  }

  if (!sequence) {
    return new Response(
      JSON.stringify({ message: "No protein sequence provided." }),
      {
        status: 400,
      }
    );
  }

  // Process the sequence
  const lines = sequence.split("\n");
  const tag = lines[0].slice(1);
  const processedSequence = lines.slice(1).join("");
  if (lines[0][0] !== ">") {
    return new Response(
      JSON.stringify({
        message: "The tag line starting with `>` is missing.",
      }),
      {
        status: 400,
      }
    );
  }

  if (!processedSequence.match(RegExp("[ACDEFGHIKLMNPQRSTVWY]+"))) {
    return new Response(
      JSON.stringify({
        message:
          "The sequence provided contains non-standard amino acids or invalid characters.",
      }),
      {
        status: 400,
      }
    );
  }

  // Here, you can process the sequence (e.g., predict methylation sites)
  console.log("Sending request to HuggingFace.");
  const startTime = Date.now();
  const resp = await fetch("https://fcharih-methylsight2.hf.space/predict", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      sequence: processedSequence,
    }),
  });
  const computeTimeInSeconds = (Date.now() - startTime) / 1000;
  console.log(
    `Computed predictions for sequence of length ${processedSequence.length} in ${computeTimeInSeconds}s.`
  );
  const { results } = await resp.json();

  return new Response(
    JSON.stringify({
      message: "Predictions successfully completed.",
      proteinName: tag,
      predictions: results,
      computeTimeInSeconds,
    }),
    {
      status: 200,
      headers: { "Content-Type": "application/json" },
    }
  );
};

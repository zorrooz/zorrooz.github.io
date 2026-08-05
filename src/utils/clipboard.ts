/**
 * Copy text to the clipboard with a hidden-textarea fallback.
 * Resolves on success; rejects when both the Clipboard API and the
 * legacy execCommand path fail (callers can surface the error).
 */
export async function copyText(text: string): Promise<void> {
  try {
    await navigator.clipboard.writeText(text)
  } catch {
    const textArea = document.createElement('textarea')
    textArea.value = text
    document.body.appendChild(textArea)
    textArea.select()
    document.execCommand('copy')
    document.body.removeChild(textArea)
  }
}

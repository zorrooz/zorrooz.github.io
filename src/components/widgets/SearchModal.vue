<template>
  <ModalOverlay @close="close">
    <div class="search-panel" role="dialog" aria-modal="true" :aria-label="t('search')">
        <div class="search-input-row d-flex align-items-center">
          <i class="fas fa-magnifying-glass search-icon"></i>
          <input
            ref="searchInput"
            v-model="keyword"
            class="search-input"
            type="text"
            inputmode="search"
            :placeholder="t('searchPlaceholder')"
          />
          <button class="search-close icon-btn" :aria-label="t('close')" @click="close">
            <i class="fas fa-times"></i>
          </button>
        </div>

        <div v-if="keyword && results.length" class="search-statusbar">
          <span class="search-count">{{ results.length }} {{ t('searchResultsLabel') }}</span>
          <span class="search-esc-hint"
            ><kbd class="search-kbd">Esc</kbd>{{ t('searchEscHint') }}</span
          >
        </div>

        <div v-if="errorMsg" class="search-status">{{ errorMsg }}</div>
        <div v-else-if="keyword && !searching && results.length === 0" class="search-status">
          <i class="fas fa-circle-info"></i> {{ t('searchNoResults') }}
        </div>

        <ul v-if="results.length" class="search-results list-unstyled mb-0">
          <li v-for="item in results" :key="item.id" class="search-result-item">
            <button class="search-result-btn w-100 text-start" @click="goToResult(item)">
              <div class="d-flex justify-content-between align-items-center gap-2">
                <span class="search-result-title">{{ item.title }}</span>
                <span class="search-result-side">
                  <span v-if="item.tags && item.tags.length" class="search-result-path">
                    {{ item.tags.slice(0, 2).join(' / ') }}
                  </span>
                  <i class="fas fa-arrow-right search-result-arrow"></i>
                </span>
              </div>
              <p class="search-result-snippet mb-0">{{ snippet(item) }}</p>
            </button>
          </li>
        </ul>
      </div>
  </ModalOverlay>
</template>

<script setup lang="ts">
import { nextTick, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRouter } from 'vue-router'
import { useSearch } from '@/composables/useSearch'
import type { SearchDoc } from '@/types'
import { toLocalePath } from '@/utils/navigation'
import { ARTICLE_ROUTE_PREFIX } from '@/config'
import ModalOverlay from '@/components/common/ModalOverlay.vue'

const emit = defineEmits(['close'])

const { t } = useI18n()
const router = useRouter()

const { keyword, results, searching, errorMsg } = useSearch()
const searchInput = ref<HTMLInputElement | null>(null)

function snippet(item: SearchDoc) {
  const text = item.description || item.content || item.title
  return text.length > 140 ? `${text.slice(0, 140)}...` : text
}

function goToResult(item: SearchDoc) {
  close()
  router.push(toLocalePath(`${ARTICLE_ROUTE_PREFIX}/${item.path}`))
}

function close() {
  emit('close')
}

nextTick(() => {
  searchInput.value?.focus()
})
</script>

<style scoped>
.search-panel {
  width: min(600px, 92vw);
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: var(--radius-lg);
  box-shadow: var(--shadow-modal);
  padding: var(--sp-4);
  max-height: 68vh;
  overflow-y: auto;
}

.search-input-row {
  gap: var(--sp-3);
  padding: var(--sp-2) var(--sp-2);
  border-bottom: 1px solid var(--line);
  padding-bottom: var(--sp-3);
}

.search-icon {
  color: var(--primary);
  font-size: 15px;
}

.search-input {
  flex: 1;
  border: none;
  outline: none;
  background: transparent;
  font-size: 16px;
  color: var(--fg);
  font-family: var(--font-sans);
}

.search-input:focus {
  box-shadow: none;
}

.search-input::placeholder {
  color: var(--fg-3);
}

.search-kbd {
  font-family: var(--font-mono);
  font-size: 11px;
  font-weight: 600;
  color: var(--fg-3);
  background: var(--surface-2);
  border: 1px solid var(--line);
  border-bottom-width: 2px;
  border-radius: var(--radius-sm);
  padding: 2px 6px;
  white-space: nowrap;
  line-height: 1.4;
}

.search-close {
  width: 30px;
  height: 30px;
  font-size: 13px;
}

.search-statusbar {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: var(--sp-3) var(--sp-2) 0;
}

.search-count {
  font-size: 11px;
  color: var(--primary);
  font-weight: 600;
  letter-spacing: 0.02em;
}

.search-esc-hint {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  font-size: 11px;
  color: var(--fg-3);
}

.search-status {
  padding: var(--sp-4) var(--sp-2) var(--sp-2);
  color: var(--fg-3);
  font-size: 14px;
  display: flex;
  align-items: center;
  gap: 8px;
}

.search-status i {
  font-size: 13px;
  color: var(--primary);
}

.search-results {
  margin-top: var(--sp-3);
}

.search-result-item + .search-result-item {
  border-top: 1px solid var(--line);
}

.search-result-btn {
  background: transparent;
  border: none;
  padding: var(--sp-3) var(--sp-2);
  border-radius: var(--radius-sm);
  color: var(--fg);
  transition: background-color 0.14s ease;
  cursor: pointer;
}

.search-result-btn:hover,
.search-result-btn:focus-visible {
  background-color: var(--surface-2);
  outline: none;
}

.search-result-title {
  font-weight: 600;
  font-size: var(--text-base);
  color: var(--fg);
  transition: color 0.14s ease;
}

.search-result-btn:hover .search-result-title,
.search-result-btn:focus-visible .search-result-title {
  color: var(--primary);
}

.search-result-path {
  font-size: 11px;
  color: var(--fg-3);
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.search-result-side {
  display: inline-flex;
  align-items: center;
  gap: var(--sp-2);
  min-width: 0;
}

.search-result-arrow {
  font-size: 11px;
  color: var(--fg-3);
  opacity: 0;
  transform: translateX(-4px);
  transition:
    opacity 0.14s ease,
    transform 0.14s ease,
    color 0.14s ease;
  flex-shrink: 0;
}

.search-result-btn:hover .search-result-arrow,
.search-result-btn:focus-visible .search-result-arrow {
  opacity: 1;
  transform: translateX(0);
  color: var(--primary);
}

.search-result-snippet {
  margin-top: 3px;
  font-size: var(--text-sm);
  color: var(--fg-2);
  display: -webkit-box;
  -webkit-line-clamp: 2;
  -webkit-box-orient: vertical;
  overflow: hidden;
  line-height: 1.5;
}
</style>
